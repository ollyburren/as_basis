library(data.table)
library(magrittr)
library(cupcake)
library(ggplot2)
library(annotSnpStats)
# aggregate chromosomes across split sample files

SHRINKAGE_METHOD<-'ws_emp'

in.dir <- '/home/ob219/rds/hpc-work/as_basis/jia_ind_analysis'
by.run <- list.files(path=in.dir,pattern='*.RDS',full.names=TRUE)
#by.run <- split(files,strsplit(basename(files),'_') %>% sapply(.,'[[',2))

## loop through compiling a list of SNPs
if(FALSE){
  snp.list<-lapply(by.run,function(x){
    tmp<-readRDS(x)
    tmp$snps$pid
  }) %>% do.call('c',.)

  ## create an alternative af file to create basis

  mtmp <- fread("/home/ob219/rds/rds-cew54-wallace-share/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab")
  td<-mtmp[pid %in% snp.list,]
  write.table(td,file="/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_jia.tab",row.names=FALSE,quote=FALSE,sep="\t")
}
## compute the basis

ref_af_file<-"/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_jia.tab"
gwas_data_dir <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)


## Clair has provided stuff based on a different set of ID's
x<-'/home/ob219/rds/hpc-work/as_basis/jia_ind_analysis/chr22.RDS'

## load in OR list


ind.proj.DT<-lapply(by.run,function(x){
  message(x)
  tmp<-readRDS(x)
  #tmp.DT<-data.table(tmp$proj.lor)[,samps,with=FALSE]
  samp.DT<-data.table(tmp$samples)
  ## get rid of unaffecteds
  samp.keep <- which(!is.na(samp.DT$alt_ilar_code))
  tmp.DT<-data.table(tmp$proj.lor[,samp.keep])
  setnames(tmp.DT,as.character(samp.DT[samp.keep,]$ID_2))
  #samp.DT[,nord:=1:.N]
  #setkey(samp.DT,ID_1)
  #setnames(tmp.DT,fam[samp.DT][order(nord),]$V2)
  ## assume that samples are in the same order but need to check !!!
  tmp.DT <- cbind(tmp.DT,pid=tmp$snps$pid)
  mDT<-melt(tmp.DT,id.vars="pid")
  setnames(mDT,c('pid','trait','beta'))
  mDT[,or:=exp(beta)]
  mDT[,.(pid,trait,or)]
}) %>% rbindlist


## unmanagable like this split into blocks of 100 samples
by.samp <- split(ind.proj.DT,ind.proj.DT$trait)
samp.groups <- split(names(by.samp),cut(1:length(names(by.samp)),50))


all.proj <- lapply(samp.groups,function(ss){
  message(length(ss))
  sample.proj.DT <- by.samp[ss] %>% rbindlist
  setkey(sample.proj.DT,pid)
  mat.emp <- create_ds_matrix(sample.proj.DT,shrink.DT,SHRINKAGE_METHOD)
  if(!identical(colnames(mat.emp),colnames(basis.mat.emp)))
    stop("Something wrong basis and projection matrix don't match")
  predict(pc.emp,newdata=mat.emp)
}) %>% do.call('rbind',.)

## assume that SNPs are in the same order for time being

control.load <- pc.emp$x["control",]

interesting.pcs <- lapply(seq_along(control.load),function(i){
  t.test(all.proj[,i],mu=control.load[i])$p.value
}) %>% do.call('c',.)

p.DT <- data.table(all.proj)
p.DT[,label:=""]
b.DT <- data.table(pc.emp$x)
b.DT[,label:=rownames(pc.emp$x)]
all.DT <- rbind(b.DT,p.DT)

library(ggplot2)
library(cowplot)
library(ggrepel)

pp<-ggplot(all.DT,aes(x=PC1,y=PC3,label=label,alpha=label!="",col=label)) + geom_point() +
geom_vline(xintercept=pc.emp$x["control",1],col='black',alpha=0.4) +
geom_hline(yintercept=pc.emp$x["control",3],col='black',alpha=0.4) +
geom_text_repel(aes(label=label)) + coord_cartesian(xlim=c(-0.05,0.05),ylim=c(-0.05,0.05)) +
scale_alpha_discrete(range=c(0.5,1)) +
theme(legend.position="none")

pp<-ggplot(all.DT,aes(x=PC1,y=PC3,label=label,alpha=label!="",col=label)) + geom_point() +
geom_vline(xintercept=pc.emp$x["control",1],col='black',alpha=0.4) +
geom_hline(yintercept=pc.emp$x["control",3],col='black',alpha=0.4) +
geom_text_repel(aes(label=label))  +
scale_alpha_discrete(range=c(0.5,1)) +
theme(legend.position="none")

save_plot("~/tmp/JIA_new.pdf",pp)
saveRDS(all.proj,file='/home/ob219/rds/hpc-work/as_basis/jia_ind_analysis/projections/jia_may.RDS')
