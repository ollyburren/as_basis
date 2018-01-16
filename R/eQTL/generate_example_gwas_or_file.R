library(data.table)

## THIS CODE CREATES A GWAS INPUT FILE

shrink.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/shrink.RDS"
pc.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/pc.RDS"
gt.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/gt_matrix.RDS"
ot.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/outcome_matrix.RDS"
snps.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/snps.RDS"
out.dir <- "/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/simulations/gwas_stats"

shrink.DT <- readRDS(shrink.file)
pc.emp <- readRDS(pc.file)
g.mat <- readRDS(gt.file)
ot.mat <- readRDS(ot.file)
snp.ids <- readRDS(snps.file)

## this is random can be any probe.idx within bounds
probe.idx<-4
bsnps <- fread("/scratch/ob219/as_basis/support_tab/as_basis_snp_support_eQTL.tab")

computeSEORMat <- function(Y,G,n0=214,n1=212){
  if(is.matrix(Y) & is.matrix(G))
    stop("Both X and Y are matrices in computeORMat")
  unx <- colSums((1-Y) * G)
  ex <- colSums(Y * G)
  a <- (n0-unx)
  b <- unx
  c <- (n1-ex)
  d <- ex
  sqrt(1/a + 1/b + 1/c + 1/d)
}

computeORMat <- function(Y,G,n0=214,n1=212){
  # ploidy is two - two groups are not identical if we use median
  if(is.matrix(Y) & is.matrix(G))
    stop("Both X and Y are matrices in computeORMat")
  unx <- colSums((1-Y) * G)
  ex <- colSums(Y * G)
  #(log(n0-unx) + log(ex)) - (log(unx) + log(n1-ex)) ##==
  log(n0-unx) + log(ex) - log(unx) - log(n1-ex)
}


NO <- data.table(pid=snp.ids,lor=computeORMat(ot.mat[,probe.idx],g.mat),se=computeSEORMat(ot.mat[,probe.idx],g.mat))[!is.infinite(lor) | lor==0,]
setkey(NO,pid)
setkey(bsnps,pid)

M <- NO[bsnps]
## for some due to the sample size we never observe one allele
## for these we fill in our own OR
M[is.na(lor),c('lor','se'):=list(0,0.44)]
M[,p.val:=pnorm(abs(lor/se),lower.tail=FALSE) * 2]
## our gwas file needs pid,a1,a2,OR,p.val
out <- M[,.(pid,ref_a1,ref_a2,or=signif(exp(lor),digits=3),p.value=signif(p.val,digits=3))]
fname <- file.path(out.dir,sprintf("probe_%s.tab",probe.idx))
write.table(out,file=fname,sep="\t",row.names=FALSE,quote=FALSE)
### code here to check our calculations
if(FALSE){
  #plot(get_qq_dt(-log10(M[!is.na(p.val),]$p.val),0.9))
  plot(get_qq_dt(-log10(M$p.val),0.9))
  abline(a=0,b=1)
  min.snp <- M[which.min(M$p.val),]
  g.mat.idx<-which(snp.ids==min.snp$pid)
  summary(glm(Y~X,dat=data.table(Y=ot.mat[,probe.idx],X=g.mat[,g.mat.idx]),family="binomial"))
  plot(head(computeSEORMat(ot.mat[,1],g.mat),n=10000),head(computeSEORMat(ot.mat[,3],g.mat),n=10000))
}
