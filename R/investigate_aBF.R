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
    exp.Z<-est_statistic(n.ctrl,n.case,snps,W=snps[CV],gamma1=CV.or,hap.freq,FP) ## NB gamma1 != gamma above
    list(exp.Z=exp.Z,CV=CV)
}

## as above but selects a constant causal variant and then estimates Z statistic based on vector of integers in sampSizeVector
## returns a list where exp.Z is a list of expectations under different sample sizes and CV is the index of the causal variant (same for all simulations)

simGWASConstantCVSampSize<-function(hap.freq,CV.or,sampSizeVector){
    snps<-head(names(hap.freq),-1)
    ## vector of odds ratios at casual variants
    CV=sample(seq_along(snps),length(CV.or))
    #message(sprintf("SELECTING CV %d",CV))
    FP <- make_GenoProbList(snps=snps,W=snps[CV],freq=hap.freq)
    tmp<-lapply(sampSizeVector,function(s){
        n.ctrl<-n.case<-s
        ## if we have exactly balanced cases and controls then approximate BF shrinkage on variance of Beta becomes infinite !
        est_statistic(n.ctrl*1.1,n.case,snps,W=snps[CV],gamma1=CV.or,hap.freq,FP) ## NB gamma1 != gamma above
    })
    names(tmp)<-make.names(sampSizeVector)
    list(exp.Z=tmp,CV=CV)
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
n.causal.blocks<-50
CV.lor<-log(1.4)
#N0 <- 3000
#N1 <- 2000
prior.abf<-10^-4
#N <- N1+N0
n.blocks<-1:length(chr1.r2s)
## sample a set of blocks
fs<-sample(n.blocks,n.causal.blocks)
## create basis of related disease
sim.name<-'sim.sample.size'
sim.base.dir<-file.path(out.dir,sim.name)
if(!dir.exists(sim.base.dir))
    dir.create(sim.base.dir)

## as regions have shared CV we need to precompute these then pick below
regions2sim<-names(chr1.r2s)[fs]

## for this disease simulate under matched cases and control numbers (optimum use of power)


## use simGWAS to compute expected z scores which we then sample from MVN
computeEZSampleSize<-function(r,CV.lor,sampleSizeVector){
    message(sprintf("Processing %s",r))
    ## for future might be worth loading all before we begin this
    ## will make things significantly quicker I should think
    load(file.path(out.dir.hap,paste0(r,'.RData')))
    simGWASConstantCVSampSize(freq,CV.lor,sampleSizeVector)
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

# posList<-chr1.r2s;sregions<-regions2sim;assoc.sim<-r2s;sigma.dir<-out.dir.sigma;maf.dir<-out.dir.maf;N0<-s*1.1;N1<-s;prior.abf=1e-4;n=1;weightByPP=TRUE;

simulateGWAS<-function(posList,sregions,assoc.sim,sigma.dir,maf.dir,N0,N1,prior.abf=1e-4,n=1,weightByPP=TRUE){
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
            pp<-apply(sim.z,1,approx.bf.z,f=maf,N=N0+N1,s=N1/(N0+N1),pi_i=prior.abf)
            lor * pp
            #pp
        }else{
            lor
        }
    })
    ## this should be in the correct order
    do.call('rbind',met)
}

## quicker to precompute expected Z's under association

## if we assume that it's optimal to have a control set 3 times case set -- need a citation for this.
#cc.size<-c(100,250,500,1000,2500,5000,10000,50000,1e5,1e6)
cc.size<-c(2e3,4e3,1e4,2e4,5e4,1e5)
simr.z<-lapply(regions2sim,computeEZSampleSize,CV.lor=CV.lor,sampleSizeVector=c(700,cc.size))
names(simr.z)<-regions2sim

## compute the basis
all.sims<-lapply(cc.size,function(s){
    message(sprintf("Processing size %s",s))
    ts<-make.names(s)
    r2s<-lapply(simr.z,function(x){
        x$exp.Z[[ts]]
    })
    simulateGWAS(chr1.r2s,regions2sim,r2s,out.dir.sigma,out.dir.maf,s*1.1,s,n=1,weightByPP=TRUE)
})



DT<-as.data.table(do.call('cbind',all.sims))
#setnames(DT,c('ss.500','ss.1e5'))
na<-make.names(cc.size)
setnames(DT,na)


## this is how these two look ploted against each other
ggplot(DT,aes(x=DT[[na[1]]],y=DT[[na[2]]])) + geom_point() + geom_abline(slope=1,color="red") + xlab(na[1]) + ylab(na[2])
##
DT$control<-0
mat<-t(DT)
pca<-prcomp(mat,scale = FALSE,center = TRUE)
tp<-melt(pca$x,id.vars=colnames(pca$x))
names(tp)<-c('size','PC','eV')

ggplot(tp,aes(x=PC,y=eV,color=size,group=size)) + geom_point() + geom_line()

## compute eucledian distance weighted by variance explained for each PC
nEucledian<-function(simLoad,actLoad,vexp){
    sqrt((actLoad-simLoad)^2 %*% vexp)
}



vexp<-summary(pca)$importance[2,]
#compute the distance between largest GWAS and smallest
comp<-rbindlist(lapply( rownames(pca$x),function(z){
    data.table(component=z,distance=nEucledian(pca$x["X1e.05",],pca$x[z,],vexp))
}))

r2s<-lapply(simr.z,function(x){
    x$exp.Z[['X1e.05']]
})
to.proj<-simulateGWAS(chr1.r2s,regions2sim,r2s,out.dir.sigma,out.dir.maf,1e5,1e5,n=100,weightByPP=FALSE)
foo<-predict(pca,t(to.proj))

bar<-apply(foo,1,nEucledian,pca$x["X1e.05",],vexp)
cbar<-apply(foo,1,nEucledian,pca$x["control",],vexp)
