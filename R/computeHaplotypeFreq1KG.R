library(devtools)
library(mvtnorm)
devtools::install_github("mdfortune/simGWAS")
library(simGWAS)



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

getNullSim<-function(sigma,n=8){
    rmvnorm(n=n,mean=rep(0,nrow(sigma)),sigma=sigma)
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

simr.z<-lapply(regions2sim,computeEZ,CV.lor=CV.lor,N0=N0,N1=N1)
names(simr.z)<-regions2sim

all.sims<-lapply(names(fs),function(d){
    message(sprintf("Processing disease %s",d))
    sdir<-file.path(sim.base.dir,d)
    if(!dir.exists(sdir))
        dir.create(sdir)
    sregions<-names(chr1.r2s)[fs[[d]]]
    nregions<-names(chr1.r2s)[!names(chr1.r2s) %in% sregions]
    toSim<-names(chr1.r2s) %in% names(chr1.r2s)[fs[[d]]]
    met<-lapply(seq_along(chr1.r2s),function(i){
        r<-names(chr1.r2s)[i]
        load(file.path(out.dir.sigma,paste0(r,'.RData')))
        load(file.path(out.dir.maf,paste0(r,'.RData')))
        if(toSim[i]){
            z<-simr.z[[r]]
        }else{
            z<-rep(0,length(maf))
        }
        sim.z<-rmvnorm(n=1,mean=z,sigma=sigma)
        lor<-estLogORfromZ(sim.z,maf,N0,N1)
        pp<-as.vector(approx.bf.z(sim.z,maf,N,N1/N0,prior.abf))
        ## this is the metric we want to return
        lor * pp
    })
    ## this should be in the correct order 
    do.call('c',met)
})



library(reshape2)
library(ggplot2)
mat<-do.call('rbind',all.sims)
mat<-rbind(mat,0)
rownames(mat)<-c(names(fs),'control')
pca<-prcomp(mat,scale = FALSE,center = TRUE)
tmp<-pca$x
mall<-melt(tmp)
names(mall)<-c('disease','pc','projection')
reo<-mall[mall$pc=='PC1',]
mall$disease<-factor(mall$disease)
## add a category so can see the projection we are comparing
mall$cat<-'simulation'
mall[mall$disease=='control',]$cat<-'control'
vexp<-summary(pca)$importance[2,]
mall$var.exp<-vexp[mall$pc]
mall$label<-sprintf("%s (%.2f)",mall$pc,mall$var.exp)
mall$label<-factor(mall$label,levels=unique(mall[order(mall$pc),]$label))
cols<-c(simulation='red',other='black')
#ggplot(mall,aes(x=label,y=projection,group=disease,color=cat)) + geom_point() + geom_path()  + theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + scale_color_manual(values=cols)

ggplot(mall,aes(x=label,y=projection,group=disease,color=disease)) + geom_point() + geom_path()  + theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + facet_grid(disease~.)

pc2.plot<-function(pc){
    all.pc<-pc$x
    DT<-data.table(all.pc[,1:2])
    DT$disease<-rownames(all.pc)
    DT$predicted=FALSE
    DT[nrow(all.pc),]$predicted<-TRUE
    ggplot(DT,aes(x=PC1,y=PC2,label=disease)) + geom_point() + geom_text(angle = 0,check_overlap=FALSE,nudge_y = 0.05) + theme_bw()
}

pc2.plot(pca)

