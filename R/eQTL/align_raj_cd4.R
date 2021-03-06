library(annotSnpStats)
library(data.table)
library(magrittr)

comp<-function(cv){
  tmp<-cv
  a<-c('A','G','C','T')
  b<-c('T','C','G','A')
  idx<-split(1:length(cv),cv)
  for(i in seq_along(a)){
    cv[idx[[a[i]]]]<-b[i]
  }
  cv
}

out_dir <- '/scratch/ob219/as_basis/raj_cd4_ind_analysis/'
data.dir <-  '/scratch/wallace/twas/'
load(file.path(data.dir,'raj-cd4-genotypes.RData'))
or.threshold <- 2
nsims <- 1e6
n.sample <- 2500
or.prior <- 0.01
support.dir <- '/scratch/ob219/as_basis/support_tab'
lor.lu.file <- sprintf("lor_posterior_%g_%g_%g_%g.tab",or.threshold,or.prior,n.sample,nsims) %>%
    file.path(support.dir,.)


out_dir <- sprintf("eqtl_ind_proj_%g_%g_%g_%g",or.threshold,or.prior,n.sample,nsims) %>%
    file.path(out_dir,.)
dir.create(out_dir, showWarnings = FALSE)

## ok add position to SNPs.

raj.snps <- snps(X) %>% data.table
raj.snps[,snp.idx:=1:nrow(raj.snps)]

library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
snps<-SNPlocs.Hsapiens.dbSNP144.GRCh37
lu<-snpsById(snps,raj.snps$snp.name,ifnotfound='drop')
lu<-data.table(as.data.frame(lu))
setkey(lu,'RefSNP_id')
setkey(raj.snps,'snp.name')
raj.snps<-raj.snps[lu][order(snp.idx),.(chromosome,snp.name,allele.1,allele.2,pos,strand,alleles_as_ambig,snp.idx)]
raj.snps<-raj.snps[!is.na(pos),pid:=paste(chromosome,pos,sep=':')]
filtered <- X[,raj.snps$snp.idx]
snps(filtered) <- as.data.frame(raj.snps)
rownames(snps(filtered)) <- snps(filtered)$snp.name


ref_af_file<-'/scratch/ob219/as_basis/support_tab/as_basis_snp_support.tab'
bsnps <- fread(ref_af_file)
setkey(bsnps,pid)
sdt <- raj.snps[,c('lid','pid'):=list(1:.N,paste(chromosome,pos,sep=':'))]
setkey(sdt,pid)
sdt <- sdt[bsnps][!(is.na(ref_a1.af) | is.na(chromosome)),]
sdt <- sdt[order(lid),]


## there might be some missing SNPs for the time being we store these as
## need to add in the to the end matrices so compatible with the basis
## or else create a new bais support file ?
new.support <- bsnps[which(bsnps$pid %in% sdt$pid),]
options(scipen=999)
write.table(new.support,file='/scratch/ob219/as_basis/support_tab/as_basis_snp_support_eQTL.tab',col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
options(scipen=0)
keep <- sdt[sdt$pid %in% bsnps$pid,]$lid
filtered <- filtered[,keep]

## save filtered annotSnpStats objects so we have for later.
saveRDS(filtered,file.path(out_dir,'raj-cd4-basis_genotypes.RDS'))

## check that thing are the same
anno.snps<-data.table(snps(filtered))
if(!identical(with(anno.snps,paste(chromosome,pos,sep=':')),sdt$pid))
stop("Something went wrong after subseting snp annotations don't line up")
## work out how we need to alter things so that we line up with reference snps
## note here we expect that outcome will be wrt to allele 2 (as snpStats)
sdt[,af.wrt.a2:=col.summary(filtered)[,"RAF"]]
setnames(sdt,c('allele.1','allele.2'),c('a1','a2'))
sdt[,flag:='unprocessed']
sdt[a1==ref_a1 & a2==ref_a2,flag:='match']
sdt[a1==ref_a2 & a2==ref_a1,flag:='flip']
sdt[a1==comp(ref_a1) & a2==comp(ref_a2),flag:='match_rc']
sdt[a1==comp(ref_a2) & a2==comp(ref_a1),flag:='flip_rc']
## Chris' code for assigning lor to individuals
dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
setnames(dat,c('f','lor.00','lor.01','lor.11'))
lor <- as.matrix(dat[,-1])  # matrix, rows indexed by 100*AF, columns by gt 00/01/11

## next get genotypes for cases and project. From snpStats 1,2,3
sm <- as(filtered,'SnpMatrix')
sm<-matrix(sm,nrow=nrow(sm),ncol=ncol(sm))
gt.lor <- lapply(seq_along(sdt$af.wrt.a2),function(j){
  if(sdt$flag[j] %in% c('flip','flip_rc')){
    g <- abs(3-as.numeric(sm[,j])) + 1
    i <- round((sdt$af.wrt.a2[j])*100)
  }else{
    g <- as.numeric(sm[,j])
    i <- round((1-sdt$af.wrt.a2[j])*100)
  }
  lor[i,g]
}) %>% do.call("rbind",.)
obj <- list(snps=sdt,proj.lor=gt.lor,samples=samples(filtered))
saveRDS(obj,file=file.path(out_dir,'raj_cd4.RDS'))
