library(data.table)

## script to gather sample counts for biobank data

data.dir <- '/home/ob219/scratch/as_basis/bb'
pheno <- fread(file.path(data.dir,'phenosummary_final_11898_18597.tsv'))
process.dir<-'/home/ob219/scratch/as_basis/gwas_stats/processed'
src.files <- list.files(path=process.dir,pattern='^bb.*RData',full.names=TRUE)

sc <- lapply(src.files,function(f){
  ## from file name get relevant data
  l<-unlist(strsplit(basename(f),':'))
  code<-gsub('\\.RData$','',l[2])
  idx<-which(pheno$Field.code==code)
  cases<-pheno[idx,]$N.cases
  controls<-pheno[idx,]$N.controls
  data.table(disease=sub('\\.RData$','',basename(f)),cases=cases,controls=controls,pmid=code)
})

sc <- rbindlist(sc)
ss<-fread("/home/ob219/scratch/as_basis/gwas_stats/sample_counts.csv")
ss<-rbind(sc,ss)
write.csv(ss,file="/home/ob219/scratch/as_basis/gwas_stats/sample_counts_bb.csv",row.names=FALSE)
