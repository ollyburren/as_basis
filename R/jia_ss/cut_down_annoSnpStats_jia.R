## need to subset annot
library(data.table)
support.dir<-'/scratch/ob219/as_basis/support_tab'
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
ref_af_file<-file.path(support.dir,'as_basis_snps_with_af.tab')
bsnps <- fread(ref_af_file)
out.dir <- '/scratch/ob219/as_basis/JIA_basis_annotSnpStats'

jfil <- list.files(path="/scratch/wallace/JIA-2017-data",pattern="annotsnpstats*",full.names=TRUE)
library(annotSnpStats)

PF <- function(fn){
  G <- get(load(fn))
  sdt <- data.table(snps(G))
  keep <- which(paste(sdt$chromosome,sdt$position,sep=':') %in% bsnps$pid)
  if(length(keep) != nrow(bsnps[chr==sdt$chromosome[1]]))
    stop("Not matching")
  G <- G[,keep]
  sdt <- data.table(snps(G))[,.(chromosome,position,allele.1,allele.2,ID)]
  samps <- data.table(samples(G))
  controls <- which(samps$phenotype==0)
  sdt[,af.wrt.a2:=col.summary(G[controls,])[,"RAF"]]
  G <- G[-controls,]
  snps(G)<-as.data.frame(sdt)
  save(G,file=file.path(out.dir,basename(fn)))
}

for(f in jfil){
  message(sprintf("Processing %s",f))
  PF(f)
}

## code for double checking things
#test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(as),sets=idx,snp.data=as(as,"SnpMatrix"))
library(data.table)
support.dir<-'/scratch/ob219/as_basis/support_tab'
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
ref_af_file<-file.path(support.dir,'as_basis_snps_with_af.tab')
bsnps <- fread(ref_af_file)
(load("/scratch/wallace/JIA-2017-data/annotsnpstats-22.RData"))
sdt <- data.table(snps(G))
keep <- which(paste(sdt$chromosome,sdt$position,sep=':') %in% bsnps$pid)
G <- G[,keep]
#test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(G),snp.data=as(G,"SnpMatrix"))
test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(G),snp.data=as(G,"SnpMatrix"))
DF<-do.call('rbind',test)
DF<- apply(DF,2,unlist)
sman<-data.table(snps(G))
sman[,c('beta','Var.beta'):=list(as.numeric(DF[,"beta"]),as.numeric(DF[,"Var.beta"]))]
controls <- which(samples(G)$phenotype==0)
sman[,af.wrt.a2:=col.summary(G[controls,])[,"RAF"]]

## Chris' code for assigning lor to individuals
lor.lu.file <- file.path(support.dir,'lor_posterior.tab')

bsnps <- fread(ref_af_file)
#bsnps <- bsnps[sample(1:nrow(bsnps),10000),]
dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
setnames(dat,c('f','lor.00','lor.01','lor.11'))
lor <- as.matrix(dat[,-1])  # matrix, rows indexed by 100*AF, columns by gt 00/01/11

## next get genotypes for cases
cases.idx<-which(samples(G)$phenotype==1)
library(magrittr)
sm <- as(G,'SnpMatrix')
sm<-matrix(sm[cases.idx,],nrow=length(cases.idx),ncol=ncol(sm))
gt.lor <- lapply(seq_along(sman$af.wrt.a2),function(j){
  i <- round((1-sman$af.wrt.a2[j])*100)
  g <- as.numeric(sm[,i])
  lor[i,g]
}) %>% do.call("rbind",.)

sman[,mean.proj.lor:=rowMeans(gt.lor)]

ggplot(sman,aes(x=beta,y=mean.proj.lor)) + geom_point() + geom_smooth()

dim(gt.lor)

mlor <- melt(gt.lor)
hist(mlor$value)
summary(rowMeans(gt.lor))
hist(rowMeans(gt.lor))
library(ggplot2)
p2 <- qplot(sman$af.wrt.a2,rowMeans(gt.lor)) + geom_smooth() + geom_hline(yintercept=0,col="red") + labs(x="AF",y="mean LOR")
p2
