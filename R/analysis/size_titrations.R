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
out_dir<-file.path('/home/ob219/scratch/as_basis/t1d_gwas_sim/input_files/size_titration/')
sim_dir<-'/home/ob219/scratch/as_basis/t1d_gwas_sim/sims/'
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
DT[,c('ss','sim','chr'):=tstrsplit(V1,'_')]
SDT<-DT[,list(N=.N),by=c('ss','sim')]

#for(i in 1:nrow(SDT)){
library(parallel)
dev.null<-mclapply(1:nrow(SDT),function(i){
  mss<-SDT[i,]$ss
  simno<-SDT[i,]$sim
  ofile<-sprintf("%s_%s.tab",mss,simno)
  message(sprintf("Processing %s",ofile))
  all<-rbindlist(lapply(file.path(sim_dir,subset(DT,ss==mss & sim==simno)$V1),readRDS))
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

SDT[,file:=sprintf("%s.tab",paste(ss,sim,sep='_'))]
SDT[,cases:=as.numeric(gsub("cases([0-9]+)control.*","\\1",ss))]
SDT[,controls:=as.numeric(gsub("cases([0-9]+)control.*","\\1",ss))]
SDT[,trait:=sprintf("T1D_%s_%s",cases,sim)]
SDT[,c('basis_trait','pmid'):=list(0,0)]
mani<-SDT[,.(trait,file,cases,controls,pmid,basis_trait)]
sim_mani_file <- '/home/ob219/scratch/as_basis/support_tab/t1d_size_titration_manifest.tab'
write.table(mani,file=sim_mani_file,row.names=FALSE,quote=FALSE)
##  next we need to create basis and do the projection
mani<-fread(sim_mani_file)



gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir)
## remove SNPs that aren't in simulations as well as T1D as this will be replaced
basis.DT<-basis.DT[(!id %in% rm_snps) & trait!='T1D',]

## loop over sim files adding - creating the basis and then projecting

## need to save a specific version of af file
new_ref_af_file<-paste('t1d_size_titration',basename(ref_af_file),sep='_')
tmp<-fread(ref_af_file)[!paste(chr,position,sep=':') %in% rm_snps,]
write.table(tmp,file=file.path(support.dir,new_ref_af_file),row.names=FALSE,quote=FALSE)

bb_traits<-fread(m_file)[grep('bb_',trait),]$trait
bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,bb_traits)[!id %in% rm_snps,]

library(parallel)
all.sims<-mclapply(1:nrow(mani),function(i){
    simname<-mani[i,]$trait
    sim.DT<-get_gwas_data(sim_mani_file,file.path(support.dir,new_ref_af_file),ld_file,out_dir,simname)
    # sometimes we get or of 0 which is impossible set these to be 1
    sim.DT[or==0 | is.infinite(or),or:=1.00001]
    sim.basis.DT<-rbind(basis.DT,sim.DT)
    setkey(sim.basis.DT,pid)
    shrink.DT<-compute_shrinkage_metrics(sim.basis.DT)
    basis.mat.emp <- create_ds_matrix(sim.basis.DT,shrink.DT,'emp')
    ## need to add control where beta is zero
    basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
    pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
    bb.mat.emp<-create_ds_matrix(bb.DT,shrink.DT,'emp')
    pred.emp <- predict(pc.emp,newdata=bb.mat.emp)
    list(loading=rbind(pc.emp$x,pred.emp),vexp=summary(pc.emp)$importance[2,])
},mc.cores=4)

names(all.sims)<-mani$trait
saveRDS(all.sims,file='/home/ob219/scratch/as_basis/t1d_gwas_sim/results/sample_size_results/19_10_17.RDS')

# next we try and plot these things

ml<-list(
  CD = 'bb_CD',
  CEL = 'bb_CEL',
  MS = 'bb_MS',
  RA = 'bb_RA',
  SLE = 'bb_SLE',

  UC = 'bb_UC'
)

g <- function(M,myml){
    M <- cbind(as.data.table(M),trait=rownames(M))
    M$compare<-"none"
    for(i in seq_along(myml)) {
        M[trait %in% c(names(myml)[i], myml[i]), compare:=names(myml)[i]]
    }
    M[trait=="control",compare:="control"]
    M
}
library(ggplot2)
pdf("~/tmp/size_titration.pdf")
lapply(SDT[order(cases,as.numeric(sim)),]$trait,function(label){
  tmpml<-ml
  tmpml[[label]] <- 'bb_T1D'
  test<-all.sims[[label]]
  fp <- g(test$loading,tmpml)
  ggplot(fp,aes(x=PC1,y=PC2,color=compare,label=trait,alpha=compare!='none')) +
  geom_point() + theme_bw() + ggtitle(label) + geom_text(show.legend=FALSE) +
  scale_alpha_discrete(guide=FALSE)
  #ggplot(fp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle(label)
})
dev.off()

## compute the var weighted distances
distances <- lapply(SDT[order(cases,as.numeric(sim)),]$trait,function(label){
  test<-all.sims[[label]]
  dist <- var_weighted_eucledian(test$vexp,test$loading[label,],test$loading['bb_T1D',])
  data.table(label=label,dist=dist)
})

distances <- rbindlist(distances)
distances[,c('disease','sample_size','iteration'):=tstrsplit(label,'_')]
distances[,sample_size:=factor(sample_size,levels=unique(sort(as.numeric(sample_size))))]
factor(distances$sample_size,levels=unique(sort(as.numeric(distances$sample_size))))

ggplot(distances,aes(x=sample_size,y=dist.V1)) + geom_boxplot() + xlab('Number of cases') + ylab('Variance weighted distance between Simulated T1D and BioBank T1D') + theme_bw()

### didn't work check that the alleles are correct wrt T1D

sim.DT<-get_gwas_data(sim_mani_file,file.path(support.dir,new_ref_af_file),ld_file,out_dir,'T1D_4000_6')
act.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,'T1D')
setnames(sim.DT,'or','sim.or')
setnames(act.DT,'or','act.or')
tp <- act.DT[sim.DT][,.(act.lor=log(act.or),sim.lor=log(sim.or))]
## definite flippage
