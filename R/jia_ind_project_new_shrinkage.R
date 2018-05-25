library(data.table)
library(magrittr)
library(cupcake)
library(ggplot2)
library(annotSnpStats)
# aggregate chromosomes across split sample files

SHRINKAGE_METHOD<-'ws_emp'

in.dir <- '/home/ob219/rds/hpc-work/as_basis/jdm_ind_analysis'
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
  write.table(td,file="/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_jdm.tab",row.names=FALSE,quote=FALSE,sep="\t")
}
## compute the basis

ref_af_file<-"/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_jdm.tab"

m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)


## Clair has provided stuff based on a different set of ID's

fam <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/jdm.fam')
setkey(fam,V1)

## load in OR list


ind.proj.DT<-lapply(by.run,function(x){
  tmp<-readRDS(x)
  #tmp.DT<-data.table(tmp$proj.lor)[,samps,with=FALSE]
  tmp.DT<-data.table(tmp$proj.lor)
  samp.DT<-data.table(tmp$samples)
  samp.DT[,nord:=1:.N]
  setkey(samp.DT,ID_1)
  setnames(tmp.DT,fam[samp.DT][order(nord),]$V2)
  ## assume that samples are in the same order but need to check !!!
  tmp.DT <- cbind(tmp.DT,pid=tmp$snps$pid)
  mDT<-melt(tmp.DT,id.vars="pid")
  setnames(mDT,c('pid','trait','beta'))
  mDT[,or:=exp(beta)]
  mDT[,.(pid,trait,or)]
}) %>% rbindlist

setkey(ind.proj.DT,pid)
## assume that SNPs are in the same order for time being
mat.emp <- create_ds_matrix(ind.proj.DT,shrink.DT,SHRINKAGE_METHOD)
if(!identical(colnames(mat.emp),colnames(basis.mat.emp)))
  stop("Something wrong basis and projection matrix don't match")
pred.emp <- predict(pc.emp,newdata=mat.emp)

control.load <- pc.emp$x["control",]

interesting.pcs <- lapply(seq_along(control.load),function(i){
  t.test(pred.emp[,i],mu=control.load[i])$p.value
}) %>% do.call('c',.)

p.DT <- data.table(pred.emp)
p.DT[,label:=""]
b.DT <- data.table(pc.emp$x)
b.DT[,label:=rownames(pc.emp$x)]
all.DT <- rbind(b.DT,p.DT)

library(ggplot2)
library(cowplot)
library(ggrepel)

pp<-ggplot(all.DT,aes(x=PC1,y=PC2,label=label,alpha=label!="",col=label)) + geom_point() +
geom_vline(xintercept=pc.emp$x["control",1],col='black',alpha=0.4) +
geom_hline(yintercept=pc.emp$x["control",2],col='black',alpha=0.4) +
geom_text_repel(aes(label=label)) + coord_cartesian(xlim=c(-0.05,0.05),ylim=c(-0.05,0.05)) +
scale_alpha_discrete(range=c(0.5,1)) +
theme(legend.position="none")

save_plot("~/tmp/JDM_new.pdf",pp)

## what happens if we regress age of onset data onto this ?

ao <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/deakin/AgeOnset.csv')

p.DT <- data.table(pred.emp)
p.DT[,ind:=rownames(pred.emp)]
setkey(ao,ID)
setkey(p.DT,ind)
# someone has an age of onset less than zero not sure what this means remove it
p.DT<-ao[p.DT][AgeOnset>0]

## fit a simple linear model and return the p.value associated with the coefficent
lapply(paste0('PC',1:11),function(pc){
  DT <- p.DT[,.(AgeOnset,pc=get(`pc`))]
  data.table(PC=pc,p.val=summary(lm(DT,formula=AgeOnset~pc))$coefficient['pc',4])
}) %>% rbindlist
