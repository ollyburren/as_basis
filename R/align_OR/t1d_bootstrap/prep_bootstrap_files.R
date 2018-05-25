## prepare bootstrap

library(magrittr)
library(snpStats)
library(data.table)
library(cupcake)

## STEP 1 consolidate all filtered and imputed files into one
## only needs to be done once for each filter. (see run_gtools.R and convert_t1d_imputed.R )

## scratch space
OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap'
SUPPORT.FILE <- file.path(OUT.DIR,'/support/as_basis_snp_support_t1t1d_bootstrap.tab')
support <- fread(SUPPORT.FILE)


if(FALSE){
  ## CODE TO CONSOLIDATE AND CHECK GENOTYPES,SUPPORT FILES ETC BEFORE BOOTSTRAP
  IN.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/tmp/snpStats'
  files <- files <- list.files(path=IN.DIR,pattern="*.RData",full.names=TRUE)
  all.obj <- lapply(files,function(x) get(load(x)))
  sm<-do.call('cbind',lapply(all.obj,'[[','sm'))
  snps<-rbindlist(lapply(all.obj,'[[','snps'))
  snps[,c('pid','lorder'):=list(paste(chr,position,sep=':'),1:.N)]
  samples<-all.obj[[1]]$samples
  rownames(sm) <- samples$sampleid
  ## add in UK regions
  ofs <- data.table(get((load("/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-gwas-data/corrected-support/WTCCC-sample-support.RData"))))
  lu<-split(ofs$b58region,ofs$affy.id)
  samples[,b58region:=unlist(lu[sampleid])]
  saveRDS(samples,file=file.path(OUT.DIR,'support/wtccc_samples.RDS'))
  ## ok next we want to align alleles with support file
  setnames(support,paste('support',names(support),sep='.'))
  setkey(support,support.pid)
  setkey(snps,pid)
  merge.snps <- support[snps]
  ## remove duplicate SNP by position
  dup.idx<-merge.snps[which(duplicated(merge.snps$support.pid)),]$lorder
  sm <- sm[,-dup.idx]
  merge.snps<-merge.snps[!lorder %in% dup.idx,]
  ## force the support file and genotypes to agree !
  by.support <- setdiff(merge.snps$support.pid,support$support.pid)
  if(length(by.support)>0){
    stop("Problem in that genotypes has more SNPs than the support - check conversion")
  }
  by.gt <- setdiff(support$support.pid,merge.snps$support.pid)
  if(length(by.gt)>0){
    #This means that after conversion (perhaps due to filtering by MAF etc) we have lost some snps the support
    #file needs updating accordingly
    support <- support[!support.pid %in% by.gt,]
  }
  if(!identical(sort(support$support.pid),sort(merge.snps$support.pid)))
    stop("Problem aligning support file with genotypes please check")
  ## next align genotypes with support - they should be aligned but check
  non.match <- merge.snps[support.ref_a1!=a1 & support.ref_a2!=a2,]
  if(nrow(non.match)>0)
    stop("GT are not aligned !")
  ## next we check to see that they are aligned by plotting MAF against each other
  merge.snps<-merge.snps[order(lorder),]
  ctrl.sum<-col.summary(sm[which(samples$t1d==1),])
  ## RAF is allele 2 freq
  plot(hexbin(1-ctrl.sum$RAF,merge.snps$support.ref_a1.af))
  ## next check samples and remove those that dont have geo data
  rm.samp.idx <- which(samples$b58region==-1)
  samples <- samples[-rm.samp.idx,]
  sm<-sm[-rm.samp.idx,]
  rownames(samples) <- samples$sampleid
  if(!identical(rownames(samples),rownames(sm)))
    stop("Samples are not aligned !")
  saveRDS(samples,'/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/wtccc_samples.RDS')
  # save updated support file
  setnames(support,gsub("support.","",names(support)))
  support<-support[,.(pid,ref_a1,ref_a2,ref_a1.af,ld.block)]
  write.table(support,file=file.path(OUT.DIR,'/support/snp_support_bootstrap_USE.tab'),sep="\t",row.names=FALSE,quote=FALSE)
  colnames(sm) <- merge.snps$support.pid
  saveRDS(sm,file=file.path(OUT.DIR,'/snpStats/t1d_wtccc_bootstrap_USE.RDS'))
  ## next create the basis and shrinkage
  snp_support_file <- file.path(OUT.DIR,'/support/snp_support_bootstrap_USE.tab')
  manifest_file <- '/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv'
  gwas_data_dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'

  basis.DT<-get_gwas_data(manifest_file,snp_support_file,gwas_data_dir,TRUE)
  ## remove t1d !
  basis.DT <- basis.DT[trait!='T1D',]
  shrink.DT<-compute_shrinkage_metrics(basis.DT)
  saveRDS(shrink.DT,file=file.path(OUT.DIR,'basis/shrink_ws.RDS'))
  ## need to add control where beta is zero
  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'ws_emp')
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  saveRDS(pc.emp,file=file.path(OUT.DIR,'basis/basis_ws.RDS'))
}

## code to generate bootstraps objects with a given configuration

library(snpStats)
library(magrittr)
library(data.table)
SIM.DIR <- '/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims_ws'
n.cases<-1900
n.controls<-1900
n.bootstraps<-500
chunk.size <- 10 # how to divide up for more efficient computation

samples <- readRDS('/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/wtccc_samples.RDS')
sm <- readRDS('/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/snpStats/t1d_wtccc_bootstrap_USE.RDS')
## chris' bootstrap function
getBoot <- function(sample.DT,status,n.boot){
  idx=which(sample.DT$t1d==status)
  regions=split(idx,sample.DT$b58region[idx])
  replicate(n.boot, sapply(regions, function(x){
    if(length(x)==1){
      return(x)
    }else{
      sample(x,replace=TRUE)
    }
  }) %>% unlist(.,use.names=FALSE))
}

## create bootstrap objects as quicker to do this first and then index
idx=c(sample(which(samples$t1d==2 & samples$b58region!=-1),n.cases),
sample(which(samples$t1d==1 & samples$b58region!=-1),n.controls))
## cutdown sm object
bs.sm <- sm[idx,]
bs.samples <- samples[idx,]
bs.samples[,cc:=t1d-1]
#each column is a bootstrap sample
bs.mat <- rbind(getBoot(bs.samples,1,n.bootstraps),getBoot(bs.samples,2,n.bootstraps))
## fit logistic
bs.samples <- as.data.frame(bs.samples)
rownames(bs.samples) <- bs.samples$sampleid
dname <- paste('wtccc',n.cases,n.controls,n.bootstraps,sep='_')
ODIR <- file.path(SIM.DIR,dname)
dir.create(ODIR,showWarnings=FALSE)
saveRDS(bs.samples,file.path(ODIR,'samples.RDS'))
saveRDS(bs.sm,file.path(ODIR,'snps.RDS'))

## next save list of bootstraps
chunks <- split(1:ncol(bs.mat),1:(n.bootstraps/chunk.size))
for(n in names(chunks)){
  mat <- bs.mat[,chunks[[n]]]
  saveRDS(mat,file=file.path(ODIR,sprintf("bs_%s.RDS",n)))
}


## adapt for smaller case controls to use just one region !
## SAMPLE PIPELINE

library(snpStats)
library(cupcake)

IN.DIR <- '/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims/wtccc_1000_1000_500'
# read in bs index
bs.idx <- readRDS(file.path(IN.DIR,'bs_48.RDS'))
samples <- readRDS(file.path(IN.DIR,'samples.RDS')) %>% as.data.frame
rownames(samples) <- samples$sampleid
snps <- readRDS(file.path(IN.DIR,'snps.RDS'))
#get pca for projections
basis <- readRDS("/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/basis.RDS")
shrink <- readRDS("/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/shrink.RDS")
support <- fread("/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/snp_support_bootstrap_USE.tab")
setkey(support,pid)
#conduct logistic fit
bs.col<-1
idx <- bs.idx[,bs.col]
bs.samples <- samples[idx,]
bs.snps <- snps[idx,]
n1 <- sum(bs.samples$cc==1)
n <- nrow(bs.samples)
rownames(bs.snps) <- rownames(bs.samples)
system.time(res <- snp.rhs.estimates(cc ~ b58region,family='binomial',data=bs.samples,snp.data=bs.snps,uncertain = TRUE))
or <- sapply(res,'[[','beta')
p <- sapply(res,function(x) 2*pnorm(abs(x$beta/sqrt(x$Var.beta)),lower.tail=FALSE))
z<-sapply(res,function(x) x$beta/sqrt(x$Var.beta))
res.dt <- data.table(pid=names(res),or=exp(or),p.value=p,trait=paste("bs",bs.col,sep='_'),n=n,n1=n1)
setkey(res.dt,pid)
## next we need to add in annotation data from support file
res.dt <- support[res.dt]
res.mat.emp<-create_ds_matrix(res.dt,shrink,'emp')
pred.emp <- predict(pc.emp,newdata=res.mat.emp)

library(cowplot)
library(ggrepel)
DT<-data.table(trait=rownames(pc.emp$x),PC1=pc.emp$x[,"PC1"],PC2=pc.emp$x[,"PC2"],basis.trait=T)
library(ggplot2)
DT.res <- data.table(trait=toupper(rownames(pred.emp)),PC1=pred.emp[,"PC1"],PC2=pred.emp[,"PC2"],basis.trait=F)
ggplot(rbind(DT,DT.res),aes(x=PC1,y=PC2,label=trait,color=basis.trait)) + geom_point() +
geom_text_repel() + scale_color_manual(guide=FALSE,values=c('firebrick','dodgerblue'))


qqnorm(z)
qqline(z)
abline(a=0,b=1,lty=2,col='red')
## s



##test code for how to do bootstrap GWAS - run on queue eventually does for 20 SNPs
test.boot <- bs.mat[,1]
sample.dat <- bs.samples[test.boot,]
snp.dat <- as(bs.sm[test.boot,],"SnpMatrix")
rownames(snp.dat) <- rownames(sample.dat)
test=snp.rhs.estimates(cc ~ b58region,family='binomial',data=sample.dat,snp.data=snp.dat,sets=1:20)
