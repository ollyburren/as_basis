##  This file takes an annotSnpStats object and removes all SNPs that are not in basis

library(data.table)

# basis reference file
ref_af_file<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/rds/user/ob219/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'

bsnps <- fread(ref_af_file)
bsnps[,c('chr','pos'):=tstrsplit(pid,':')]
out.dir <- '/rds/user/ob219/hpc-work/as_basis/JIA_basis_annotSnpStats'

jfil <- list.files(path="/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/JIA-2017-data/",pattern="*.RData",full.names=TRUE)
library(annotSnpStats)

PF <- function(fn){
  G <- get(load(fn))
  #chr <- gsub("annotsnpstats-([^.]+)\\.RData","\\1",basename(fn))
  sdt <- data.table(snps(G))
  keep <- which(paste(sdt$chromosome,sdt$position,sep=':') %in% bsnps$pid)
  missing <- nrow(bsnps[chr==sdt$chromosome[1]]) - length(keep)
  if(missing!=0)
    message(sprintf("Not matching missing %d",missing))
  G <- G[,keep]
  sdt <- data.table(snps(G))[,.(position,a0=allele.1,a1=allele.2,rs_id=ID,chromosome)]
  samps <- data.table(samples(G))
  ## we don't have any control data so cannot compute RAF for them
  #controls <- which(samps$phenotype==0)
  sdt[,af.wrt.a2:=col.summary(G)[,"RAF"]]
  #G <- G[-controls,]
  snps(G)<-as.data.frame(sdt)
  alleles(G)<-c('a0','a1')
  save(G,file=file.path(out.dir,basename(fn)))
}

for(f in jfil){
  message(sprintf("Processing %s",f))
  PF(f)
}
