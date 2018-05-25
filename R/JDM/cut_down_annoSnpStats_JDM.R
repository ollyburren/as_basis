## need to subset annot
library(data.table)

support.dir<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
# reference allele frequencies
ref_af_file<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
#ld_file<-file.path(support.dir,'all.1cM.tab')
#m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv'
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/rds/user/ob219/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'

#support.dir<-'/scratch/ob219/as_basis/support_tab'
#gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
#ref_af_file<-file.path(support.dir,'as_basis_snps_with_af.tab')
bsnps <- fread(ref_af_file)
bsnps[,c('chr','pos'):=tstrsplit(pid,':')]
out.dir <- '/rds/user/ob219/hpc-work/as_basis/JDM_basis_annotSnpStats'

jfil <- list.files(path="/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/imputed-1kg-apr2018",pattern="*.RData",full.names=TRUE)
library(annotSnpStats)

PF <- function(fn){
  G <- get(load(fn))
  chr <- gsub("chr\\-([^.]+)\\.RData","\\1",basename(fn))
  sdt <- data.table(snps(G))[,chromosome:=chr]
  keep <- which(paste(sdt$chromosome,sdt$position,sep=':') %in% bsnps$pid)
  missing <- nrow(bsnps[chr==sdt$chromosome[1]]) - length(keep)
  if(missing!=0)
    message(sprintf("Not matching missing %d",missing))
  G <- G[,keep]
  sdt <- data.table(snps(G))[,.(position,a0,a1,rs_id)][,chromosome:=chr]
  samps <- data.table(samples(G))
  ## we don't have any control data so cannot compute RAF for them
  #controls <- which(samps$phenotype==0)
  sdt[,af.wrt.a2:=col.summary(G)[,"RAF"]]
  #G <- G[-controls,]
  snps(G)<-as.data.frame(sdt)
  save(G,file=file.path(out.dir,basename(fn)))
}

for(f in jfil){
  message(sprintf("Processing %s",f))
  PF(f)
}

## code for double checking things
#test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(as),sets=idx,snp.data=as(as,"SnpMatrix"))


ggplot(sman,aes(x=beta,y=mean.proj.lor)) + geom_point() + geom_smooth()

dim(gt.lor)

mlor <- melt(gt.lor)
hist(mlor$value)
summary(rowMeans(gt.lor))
hist(rowMeans(gt.lor))
library(ggplot2)
p2 <- qplot(sman$af.wrt.a2,rowMeans(gt.lor)) + geom_smooth() + geom_hline(yintercept=0,col="red") + labs(x="AF",y="mean LOR")
p2
