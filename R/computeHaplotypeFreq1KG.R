library(devtools)
library(mvtnorm)
devtools::install_github("mdfortune/simGWAS")
library(simGWAS)
library(data.table)



## tabix bin

getRegionVCF<-function(region,vcf.file,tabix_bin,pos.keep){
    chrom <- gsub("^([^:]+).*","\\1",region)
    ## get the header
    my.pipe<-pipe(paste(tabix_bin,'-H',vcf.file,region))
    header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
    close(my.pipe)
    cnames<-unlist(strsplit(header,"\t"))
    tmp<-as.data.frame(fread(paste(tabix_bin,vcf.file,region,' | grep PASS'),sep="\t",header=FALSE,stringsAsFactors=FALSE))
    colnames(tmp)<-cnames
    idx<-which(tmp$POS %in% pos.keep)
    if(length(idx) != length(pos.keep)){
        cat(sprintf("Warning %d pos.keep SNPs(%d) not found\n",length(pos.keep),length(pos.keep)-length(idx)))
    }
    tmp<-tmp[idx,]
    tmp<-tmp[order(tmp$POS),]
    gt<-tmp[,10:ncol(tmp)]
}

getSM<-function(gt){
    ## convert gt to snpMatrix object
    sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
    sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
    sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
    ## set anything else to a missing value
    sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
    sm<-new("SnpMatrix", sm)
    rownames(sm)<-1:nrow(sm)
    colnames(sm)<-1:ncol(sm)
    return(sm)
}

getHapFreq<-function(gt){
    ## for the simulation we only care about one haplotype per individual
    hap<-function(x){
        l<-strsplit(x,"")
        as.numeric(do.call('c',lapply(l,'[[',1)))
    }
    haps<-t(apply(gt,2,hap))
    rownames(haps)<-NULL
    snps <- colnames(haps) <- paste0("s",1:ncol(haps))
    ## simGWAS expects SNPHAP format
    freq <- as.data.frame(haps+1)
    hapstr<-apply(freq,1,paste0,collapse="")
    freq$hapstr<-hapstr
    freq<-freq[!duplicated(freq$hapstr),]
    thapstr<-table(hapstr)
    freq$hapcount<-thapstr[match(freq$hapstr,names(thapstr))]
    ## total haptypes is same as individuals ncols
    freq$Probability<-freq$hapcount/ncol(gt)
    freq$hapstr<-NULL
    freq$hapcount<-NULL
    freq
}

## compute covariance matrix based on LD structure
getSigma<-function(sm){
    LD<-snpStats::ld(sm,sm,stat="R",symmetric=TRUE)
    LD[which(is.na(LD))]<-0
    make.positive.definite(LD)
}


## hap.freq - use getHapFreq to generate haplotype freq
## CV.or - vector of effect sizes for causal variant
## n.ctrl - number of controls
## n.case - number of cases
simGWASr<-function(hap.freq,CV.or,n.ctrl,n.case){
    snps<-head(names(hap.freq),-1)
    N <- n.ctrl+n.case
    ## vector of odds ratios at casual variants
    CV=sample(seq_along(snps),length(CV.or))
    #message(sprintf("SELECTING CV %d",CV))
    FP <- make_GenoProbList(snps=snps,W=snps[CV],freq=hap.freq)
    est_statistic(n.ctrl,n.case,snps,W=snps[CV],gamma1=CV.or,hap.freq,FP) ## NB gamma1 != gamma above
}



## this code works out which snps are in a given LD block. It only needs to run once

if(!file.exists("/Users/oliver/DATA/AS_BASIS/ld.block.2.basis.snp.RData")){
    library(data.table)
    library(GenomicRanges)
    (load("/Users/oliver/DATA/AS_BASIS/all_or_shared_with_af.RData"))
    setkey(final.t,id)
    s.gr<-with(unique(final.t),GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,end=position)))
    ## load in LD blocks and work out overlap
    ld<-fread("/Users/oliver/DATA/AS_BASIS/all.1cM.bed")
    ld$start<-as.numeric(gsub("([^\\-]+)\\-.*","\\1",ld$V2))
    ld$end<-as.numeric(gsub(".*\\-(.*)","\\1",ld$V2))
    ld.gr<-with(ld,GRanges(seqnames=Rle(V1),ranges=IRanges(start=start,end=end)))
    ol<-as.matrix(findOverlaps(ld.gr,s.gr))
    r2s<-data.table(region=paste(ld$V1[ol[,1]],ld$V2[ol[,1]],sep=':'),position=start(s.gr)[ol[,2]])
    chr1.r2s<-r2s[grep("^1:",r2s$region),]
    chr1.r2s<-split(chr1.r2s$position,chr1.r2s$region)
    chr1.r2s<-chr1.r2s[sapply(chr1.r2s,length)>=10]
    save(chr1.r2s,file="/Users/oliver/DATA/AS_BASIS/ld.block.2.basis.snp.RData")
}else{
    load("/Users/oliver/DATA/AS_BASIS/ld.block.2.basis.snp.RData")
}


## this only needs to be run once to provide a set of LD blocks + association that can be used over all downstream
## simulations
sim.name<-'base_simulation'
out.dir<-file.path('/Users/oliver/DATA/AS_BASIS/',sim.name)
out.dir.hap<-file.path(out.dir,'chr1_hap_prob')
out.dir.sigma<-file.path(out.dir,'chr1_null_sigma')
out.dir.maf<-file.path(out.dir,'chr1_maf')
if(!dir.exists(out.dir)){
    library(data.table)
    
    for(cdir in c(out.dir,out.dir.sigma,out.dir.hap, out.dir.maf)){
        if(!dir.exists(cdir))
           dir.create(cdir)
    }
    kg.vcf.file<-'/Users/oliver/DATA/1KG/VCF/EUR/by.chr.phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz'
    tabix.bin<-'/usr/local/bin/tabix'
    bad.locus<-integer()
    for(i in seq_along(chr1.r2s)){
        pos<-chr1.r2s[[i]]
        t.region<-names(chr1.r2s)[i]
        if(length(pos)<10){
            message(sprintf("Odd this region has %d SNPs which is below threshold. Skipping",length(pos)))
            bad.locus<-c(bad.locus,i)
            next
        }
        message(paste("Processing",t.region))
        gt<-getRegionVCF(t.region,kg.vcf.file,tabix.bin,pos)
        ## compute haplotype frequencies for simulating under alternatives
        freq<-getHapFreq(gt)
        o.f<-file.path(out.dir.hap,paste0(t.region,'.RData'))
        save(freq,file=o.f)
        ## simulate under the null of no association do 100 for each region
        sm<-getSM(gt)
        sigma<-getSigma(sm)
        o.f<-file.path(out.dir.sigma,paste0(t.region,'.RData'))
        save(sigma,file=o.f)
        #null.sims<-getNullSim(sm,100)
        #o.f<-file.path(out.dir.null,paste0(t.region,'.RData'))
        #save(null.sims,file=o.f)
        ## save vector of MAFs 
        maf<-col.summary(sm)[['MAF']]
        o.f<-file.path(out.dir.maf,paste0(t.region,'.RData'))
        save(maf,file=o.f)
    }
}
## work out which are valid LD blocks (require >10 SNPs per LD Block). Want to simulate
## current basis so use same SNPs do for chromosome 1

## chromosome 1

## select groups of LD blocks such that there are 8 groups (representing a disease)


## contains code for computing overlap sets work out how to do properly 
## with environment
source("~/gitr/as_basis/R/sim_sets.R")
## this contains code for computing Wakefield's approx BF and thus posterior probabilities
source("~/gitr/as_basis/R/wakefield.R")
out.dir.hap<-file.path(out.dir,'chr1_hap_prob')
out.dir.null<-file.path(out.dir,'chr1_null_sim')
out.dir.maf<-file.path(out.dir,'chr1_maf')
n.diseases<-8
n.causal.blocks<-50
CV.lor<-log(1.4)
N0 <- 3000
N1 <- 2000
prior.abf<-10^-4
N <- N1+N0
n.blocks<-1:length(chr1.r2s)
## create basis of related disease 
fs<-selectGroups(n.blocks,n.diseases,n.causal.blocks)
names(fs)<-paste0('d',1:n.diseases)
sim.name<-'sim.1'
sim.base.dir<-file.path(out.dir,sim.name)
if(!dir.exists(sim.base.dir))
    dir.create(sim.base.dir)

## as regions have shared CV we need to precompute these then pick below
regions2sim<-names(chr1.r2s)[unique(do.call('c',fs))]

## use simGWAS to compute expected z scores which we then sample from MVN
computeEZ<-function(r,CV.lor,N0,N1){
    message(sprintf("Processing %s",r))
    ## for future might be worth loading all before we begin this 
    ## will make things significantly quicker I should think
    load(file.path(out.dir.hap,paste0(r,'.RData')))
    simGWASr(freq,CV.lor,N0,N1)
}

## think this is correct need to double check with Chris
estLogORfromZ<-function(z,maf,N0,N1){
    mt<-z/sqrt((2*maf*(1-maf)))
    as.vector(mt * sqrt(1/N0 + 1/N1))
}

# posList list with names representing a region and elements positions of SNPs within. 
# sregions - list of regions in the above to simulate association
# sigma.dir - location of precomputed covariance matrix that maps to posList 
# maf.dir - location of precompute MAF that maps to posList
# N0 - Number of cases
# N1 - Number of controls

simulateGWAS<-function(posList,sregions,assoc.sim,sigma.dir,maf.dir,N0,N1,prior.abf=10^-4,n=1,weightByPP=TRUE){
    #sregions<-names(chr1.r2s)[fs[[d]]]
    #nregions<-names(chr1.r2s)[!names(chr1.r2s) %in% sregions]
    toSim<-names(posList) %in% sregions
    met<-lapply(seq_along(posList),function(i){
        r<-names(posList)[i]
        ssfile<-file.path(sigma.dir,paste0(r,'.RData'))
        mafile<-file.path(maf.dir,paste0(r,'.RData'))
        #message(sprintf("%s,%s,%s",ssfile,mafile,r))
        sigma<-get(load(ssfile))
        maf<-get(load(mafile))
        ## compute expected Z scores
        if(toSim[i]){
            e.z<-assoc.sim[[r]]
        }else{
            e.z<-rep(0,length(maf))
        }
        sim.z<-rmvnorm(n=n,mean=e.z,sigma=sigma)
        ## note that sim.z will now be a matrix
        lor<-apply(sim.z,1,estLogORfromZ,maf=maf,N0=N0,N1=N1)
        #lor<-estLogORfromZ(sim.z,maf,N0,N1)
        #pp<-as.vector(approx.bf.z(sim.z,maf,N,N1/N0,prior.abf))
        if(weightByPP){
            pp<-apply(sim.z,1,approx.bf.z,f=maf,N=N0+N1,s=N1/N0,pi_i=prior.abf)
            lor * pp
        }else{
            lor
        }
    })
    ## this should be in the correct order 
    do.call('rbind',met)
}

## quicker to precompute expected Z's under association
simr.z<-lapply(regions2sim,computeEZ,CV.lor=CV.lor,N0=N0,N1=N1)
names(simr.z)<-regions2sim

## compute the basis
all.sims<-lapply(names(fs),function(d){
    message(sprintf("Processing disease %s",d))
    reg2sim<-names(chr1.r2s)[fs[[d]]]
    simulateGWAS(chr1.r2s,reg2sim,simr.z,out.dir.sigma,out.dir.maf,N0,N1,n=1)
})
mat<-do.call('rbind',lapply(all.sims,as.vector))
mat<-rbind(mat,0)
rownames(mat)<-c(names(fs),'control')
pca<-prcomp(mat,scale = FALSE,center = TRUE)



## given a disease compute a set of simulated gwas say 1000 times and we can use these to estimate 
sim.d<-sample(names(fs),1)
reg2sim<-names(chr1.r2s)[fs[[sim.d]]]
to.project<-simulateGWAS(chr1.r2s,reg2sim,simr.z,out.dir.sigma,out.dir.maf,N0,N1,n=1000,weightByPP = FALSE)
project.pc<-predict(pca,t(to.project))
## long and thin
fdt<-data.table(project.pc)
melt(fdt,id.vars=colnames(fdt))
fdt<-melt(fdt,measure.vars=colnames(fdt))
act<-melt(pca$x[sim.d,],measure.vars=colnames(pca$x))
act<-melt(pca$x,measure.vars=colnames(pca$x))
act$actual<-act$Var1 %in% sim.d

vexp<-summary(pca)$importance[2,]
fdt$variable<-sprintf("%s (%.2f)",fdt$variable,vexp[fdt$variable])
act$Var2<-sprintf("%s (%.2f)",act$Var2,vexp[act$Var2])


#act$variable<-rownames(act)
#ggplot(fdt,aes(x=variable,y=value)) + geom_boxplot() + geom_point(data=act,aes(x=Var2,color=actual),pch=16,size=2) + geom_text(data=act,aes(x=Var2,color=actual,label=Var1),nudge_x = 0.2) + theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

## computing a p dimensional distance.

nEucledian<-function(simLoad,actLoad){
    sqrt(sum(((actLoad-simLoad)^2)))
}

## compute single distance stat for each disease to see if we can tell what is being simulated
comp<-rbindlist(lapply( rownames(pca$x),function(z){
    data.table(component=z,distance=apply(project.pc,1,nEucledian,pca$x[z,]))
}))

ylabel=sprintf("Euc. Dist. of %s simulated data\nfrom actual disease loading",sim.d)

ggplot(comp,aes(x=component,y=distance)) + geom_boxplot() + theme_bw() + ylab("Eucledian Distance from  by PC Loading") + xlab("Disease") + ylab(ylabel)

## titrate the # initially keep the number of cases the same

#n.cases<-c(100,200,350,500,800,1000,1500,2000)
n.cases<-seq(100,2000,100)
n.sims<-100
titr<-lapply(n.cases,function(n.cas){
    message(sprintf("Simulating %d cases",n.cas))
    to.project<-simulateGWAS(chr1.r2s,reg2sim,simr.z,out.dir.sigma,out.dir.maf,N0,n.cas,n=n.sims,weightByPP = FALSE)
    project.pc<-predict(pca,t(to.project))
    rbindlist(lapply( rownames(pca$x),function(z){
        data.table(n.cases=n.cas,disease=z,distance=apply(project.pc,1,nEucledian,pca$x[z,]))
    }))
})

titr<-rbindlist(titr)

#ggplot(foo,aes(x=factor(n.cases),y=distance)) + geom_boxplot() + theme_bw() + ylab("Eucledian Distance from  by PC Loading") + xlab("Number of Cases") + ylab(ylabel) + facet_grid(component~.)
sum<-titr[,list(mdist=mean(distance),sd=sd(distance)),by=c('n.cases','disease')]
sum$me<-(sum$sd/sqrt(n.sims))*1.96
limits <- aes(ymax = mdist+me,ymin=mdist-me)

ggplot(sum,aes(x=n.cases,y=mdist,group=disease,colour=disease)) + geom_point() + geom_line() + geom_errorbar(limits, width=0.2) + theme_bw()  + ylab(sprintf("Mean Eucledian Distance from simulation(%s x %d)\nto disease  by PC Loading",sim.d, n.sims)) + xlab("Number of Cases")


## randomly sample GWAS 
## compute the expected z scores under the alternative using unimodal effect size
simr.rnd<-lapply(names(chr1.r2s),computeEZ,CV.lor=CV.lor,N0=N0,N1=N1)
names(simr.rnd)<-names(chr1.r2s)

## simulate 100 GWAS to start
n.sample.gwas<-500

## randomly select 50 regions to be causal (n.causal.blocks)
rand2sim<-sample(names(chr1.r2s),n.causal.blocks)
rand.gwas.to.project<-lapply(1:n.sample.gwas,function(i){
    message(sprintf("Processing %d",i))
    simulateGWAS(chr1.r2s,rand2sim,simr.rnd,out.dir.sigma,out.dir.maf,N0,N1,n=1)
})

foo<-do.call('cbind',rand.gwas.to.project)
foo.pc<-predict(pca,t(foo))

foo.euc<-rbindlist(lapply( rownames(pca$x),function(z){
    data.table(disease=z,distance=apply(foo.pc,1,nEucledian,pca$x[z,]))
}))


ggplot(foo.euc,aes(x=disease,y=distance)) + geom_boxplot() + theme_bw()

## w weights for a given pc and basis (eigenvector)

## principal 
w<-pca$rotation[,1]
d1<-mat['d1',]
d2<-foo[,1]
# e.distance.slow<-function(w,d1,d2){
#     DT<-data.table(w=w,d1=d1,d2=d2)
#     DT[,`:=`(w.sq=w^2),]
#     summ.DT<-DT[,list(t1=sum(w.sq*(d1^2)),t2=sum(w.sq*(d2^2)),t3=-2*sum(w.sq*d1*d2)),]
#     t4<-sum(do.call('c',lapply(1:nrow(DT),function(i){
#         tmpDT<-DT[-i,]
#         oDT<-DT[i,]
#         sum(tmpDT$w * tmpDT$d1 * oDT$w * oDT$d1)
#     })))
#     
#     t5<-sum(do.call('c',lapply(1:nrow(DT),function(i){
#         tmpDT<-DT[-i,]
#         oDT<-DT[i,]
#         sum(tmpDT$w * tmpDT$d2 * oDT$w * oDT$d2)
#     })))
#     
#     t6<--2*sum(do.call('c',lapply(1:nrow(DT),function(i){
#         tmpDT<-DT[-i,]
#         oDT<-DT[i,]
#         sum(tmpDT$w * tmpDT$d2 * oDT$w * oDT$d1)
#     })))
#     sum(summ.DT$t1,summ.DT$t2,summ.DT$t3,t4,t5,t6)
# }
# 
# e.distance.fast1<-function(w,d1,d2){
#     DT<-data.table(w=w,d1=d1,d2=d2)
#     mat<-matrix(0,nrow=length(w),ncol=6)
#     mat[,1]<-DT$w^2 * DT$d1^2
#     mat[,2]<-DT$w^2 * DT$d2^2
#     mat[,3]<-DT$w^2 *  DT$d1 *  DT$d2
#     for(i in 1:length(w)){
#         tmpDT<-DT[-i,]
#         oDT<-DT[i,]
#         mat[i,4]<-sum(tmpDT$w * tmpDT$d1 * oDT$w * oDT$d1)
#         mat[i,5]<-sum(tmpDT$w * tmpDT$d2 * oDT$w * oDT$d2)
#         mat[i,6]<-sum(tmpDT$w * tmpDT$d2 * oDT$w * oDT$d1)
#     }
#     sum(colSums(mat) * c(1,1,-2,1,1,-2))
# }
# 
# e.distance.fast2<-function(w,d1,d2){
#     dat<-do.call('cbind',list(w,d1,d2))
#     mat<-matrix(0,nrow=length(w),ncol=6)
#     mat[,1]<-dat[,1]^2 * dat[,2]^2
#     mat[,2]<-dat[,1]^2 * dat[,3]^2
#     mat[,3]<-dat[,1]^2 * dat[,2] * dat[,3]
#     for(i in 1:length(w)){
#         ## get line we are interested in
#         odat<-dat[i,]
#         dat[i,]<-0
#         #mat[i,4]<-sum(dat[,1] * dat[,2] * odat[1] * odat[2])
#         dat[,1]
#         mat[i,5]<-sum(dat[,1] * dat[,3] * odat[1] * odat[3])
#         mat[i,6]<-sum(dat[,1] * dat[,3] * odat[1] * odat[2])
#         dat[i,]<-odat
#     }
#     sum(colSums(mat) * c(1,1,-2,1,1,-2))
# }

library(Rcpp)

cppFunction('NumericVector computeTerms(NumericMatrix x){
    int nrow = x.nrow();
    NumericVector out(3);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < nrow; j++) {
            if(i==j) continue;
            out[0] += x(j,0) * x(j,1) * x(i,0) * x(i,1);
            out[1] += x(j,0) * x(j,2) * x(i,0) * x(i,2);
            out[2] += x(j,0) * x(j,2) * x(i,0) * x(i,1);
        }
    }
    return(out);
}')

e.distance.fast3<-function(w,d1,d2){
    dat<-do.call('cbind',list(w,d1,d2))
    mat<-matrix(0,nrow=length(w),ncol=3)
    mat[,1]<-dat[,1]^2 * dat[,2]^2
    mat[,2]<-dat[,1]^2 * dat[,3]^2
    mat[,3]<-dat[,1]^2 * dat[,2] * dat[,3]
    sum(c(colSums(mat),computeTerms(dat)) * c(1,1,-2,1,1,-2))
}

## is this the quickest way to do things ?


# n<-1000
# library(microbenchmark)
# microbenchmark(
#     e.distance.fast1(w[1:n],d1[1:n],d2[1:n]),
#     e.distance.fast2(w[1:n],d1[1:n],d2[1:n]),
#     e.distance.fast3(w[1:n],d1[1:n],d2[1:n])
# )

## these should be equivalent
(sum(DT$w*DT$d1) - sum(DT$w*DT$d2))^2
e.distance.fast3(DT$w,DT$d1,DT$d2)
## now need to do it across all pc for a given disease
sqrt(sum(apply(pca$rotation,2,e.distance.fast3,d1=d1,d2=d2)))



## theoretical distance calculations (for maths see Overleaf )



# MahalanobisDistance<-function(A,B){
#     n1 <- length(A)
#     n2 <- length(B)
#     n <- n1 + n2
#     xdiff<-mean(A)-mean(B)
#     cA <- cov(A)
#     cB <- cov(B)
#     pC <- n1/n*cA + n2/n^cB
#     sqrt(xdiff * solve(pC) * t(xdiff))
# }
# 
# tmA<-matrix(c( 2, 2, 2, 5, 6, 5, 7, 3, 4, 7, 6, 4, 5, 3, 4, 6, 2, 5, 1, 3),ncol=2,byrow=TRUE)
# tmB<-matrix(c( 6, 5, 7, 4, 8, 7, 5, 6, 5, 4),ncol=2,byrow=TRUE)
# 
# MahalanobisDistance(tmA,tmB)


# pc2.plot<-function(pc,...){
#     all.pc<-pc$x
#     DT<-data.table(all.pc[,1:2])
#     DT$disease<-rownames(all.pc)
#     DT$predicted=FALSE
#     args <- list(...)
#     if("proj" %in% names(args)){
#         pDT<-data.table(args[['proj']][,1:2])
#         #pDT$disease<-paste0('p',1:nrow(pDT))
#         pDT$disease<-''
#         pDT$predicted=TRUE
#         DT<-rbind(DT,pDT)
#     }
#     ggplot(DT,aes(x=PC1,y=PC2,label=disease,colour=predicted)) + geom_point() + geom_text(angle = 0,check_overlap=FALSE,nudge_y = 0.05) + theme_bw()
# }
# 
# computeDist<-function(pc,...){
#     all.pc<-pc$x
#     DT<-all.pc[,1:2]
#     args <- list(...)
#     if("proj" %in% names(args)){
#         pDT<-args[['proj']][,1:2]
#         #pDT$disease<-paste0('p',1:nrow(pDT))
#         DT<-rbind(DT,pDT)
#     }
#     dist(DT)
# }
# 
# pc2.plot(pca,proj=project.pc)
# 
# pdist<-computeDist(pca,proj=project.pc)
# 
