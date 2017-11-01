library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)

## import GWAS data for basis
## support files
support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snps_with_alleles.tab')
ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-file.path(support.dir,'as_basis_manifest.tab')
out_dir<-file.path('/home/ob219/scratch/as_basis/t1d_gwas_sim/input_files/null/')
sim_dir<-'/home/ob219/scratch/as_basis/t1d_gwas_sim/null_sims/'
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'

## generate this as follows
#foo<-fread("/home/ob219/scratch/as_basis/t1d_gwas_sim/input_files/size_titration/cases800control800_7.tab")
#foo$id<-paste(foo$chr,foo$position,sep=':')
#bas<-fread('/home/ob219/scratch/as_basis/support_tab/as_basis_snps.tab')
#bas$id<-paste(bas$chr,bas$position,sep=':')
#idx<-which(bas$id %in% foo$id)
#write(bas[-idx,]$id,file='/home/ob219/scratch/as_basis/support_tab/rm_snps_size_titration.txt')


snp_exlude_file<-'/home/ob219/scratch/as_basis/support_tab/rm_snps_size_titration.txt'
rm_snps<-scan(snp_exlude_file,character())
bsnps<-fread(ref_af_file)
bsnps[,pid:=paste(chr,position,sep=':')]
setkey(bsnps,pid)
## TEST

DT<-data.table(list.files(path=sim_dir,pattern='*.RDS'))
DT[,c('disease','ss','sim','chr'):=tstrsplit(V1,'_')]
SDT<-DT[,list(N=.N),by=c('sim')]

#for(i in 1:nrow(SDT)){
library(parallel)
dev.null<-mclapply(1:nrow(SDT),function(i){
  simno<-SDT[i,]$sim
  ofile<-sprintf("%s.tab",simno)
  message(sprintf("Processing %s",ofile))
  all<-rbindlist(lapply(file.path(sim_dir,subset(DT,sim==simno)$V1),readRDS))
  all[,pid:=paste(chr,position,sep=':')]
  ## need to check that all basis SNPs are available
  setkey(all,pid)
  res<-all[bsnps]
  if(nrow(res)!=nrow(bsnps)){
    stop(sprintf("Number of SNPs are missing got %d had %d",nrow(res),nrow(bsnps)))
  }
  not.found<-res[is.na(id),]$pid
  res<-res[order(chr,position),.(id,chr,position,p.val,or,a1,a2)][!is.na(id),]
  res<-res[is.na(or),c('p.val','or'):=list(0.9999,1.0000001)]
  res<-align_alleles(res,bsnps[!pid %in% not.found,],check=FALSE)
  ## or are flipped thus we need 1/or !!
  res<-res[,.(id,chr,position,p.val=sprintf("%.4f",p.val),or=sprintf("%.4f",1/or))]
  options(scipen=999)
  write.table(res,file=file.path(out_dir,ofile),sep="\t",row.names=FALSE,quote=FALSE)
  options(scipen=0)
},mc.cores=4)

## create a manifest file so we can load in using cupcake

SDT[,file:=sprintf("%s.tab",sim)]
SDT[,cases:=5913]
SDT[,controls:=5919]
SDT[,trait:=sprintf("T1D_%s",sim)]
SDT[,c('basis_trait','pmid'):=list(0,0)]
mani<-SDT[,.(trait,file,cases,controls,pmid,basis_trait)]
sim_mani_file <- '/home/ob219/scratch/as_basis/support_tab/null_manifest.tab'
write.table(mani,file=sim_mani_file,row.names=FALSE,quote=FALSE)
##  next we need to create basis and do the projection
mani<-fread(sim_mani_file)



gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir)
## remove SNPs that aren't in simulations as well as T1D as this will be replaced
basis.DT<-basis.DT[(!id %in% rm_snps) ,]
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
## loop over sim files adding - creating the basis and then projecting

## need to save a specific version of af file
new_ref_af_file<-paste('null',basename(ref_af_file),sep='_')
tmp<-fread(ref_af_file)[!paste(chr,position,sep=':') %in% rm_snps,]
write.table(tmp,file=file.path(support.dir,new_ref_af_file),row.names=FALSE,quote=FALSE)

bb_traits<-fread(m_file)[grep('bb_',trait),]$trait
bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,bb_traits)[!id %in% rm_snps,]

## this is different need to project the null onto the the actual basis to see how far null phenotype sits


library(parallel)
all.sims<-mclapply(1:nrow(mani),function(i){
    simname<-mani[i,]$trait
    sim.DT<-get_gwas_data(sim_mani_file,file.path(support.dir,new_ref_af_file),ld_file,out_dir,simname)
    # sometimes we get or of 0 which is impossible set these to be 1
    sim.DT[or==0 | is.infinite(or),or:=1.00001]
    sim.mat.emp<-create_ds_matrix(sim.DT,shrink.DT,'emp')
    data.table(predict(pc.emp,newdata=sim.mat.emp))
},mc.cores=4)

all.sims<-rbindlist(all.sims)

act<-pc.emp$x
tmp.rn<-rownames(act)
act<-data.table(act)
all<-rbind(act[,c('trait','sim'):=list(tmp.rn,FALSE)],
  all.sims[,c('trait','sim'):=list(mani$trait,TRUE)])
plot <- melt(all,id.vars=c('trait','sim'))
## make control a sim
#plot[trait=='control',sim:=TRUE]
plot[,label:=ifelse(sim,'simulation',trait)]

## plot scree plots
library(ggplot2)
ggplot(plot,aes(x=variable,y=value,alpha=!sim,color=label,group=trait)) + geom_point() + geom_line() + theme_bw() + scale_alpha_discrete(range = c(0.5, 1))
