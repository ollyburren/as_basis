## next we compute GWAS at different sample size configurations and store the
## summary statistics

library(snpStats)
library(data.table)
library(cupcake)

SUPPORT.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/snp_support_bootstrap_USE.tab"

samples <- readRDS('/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/wtccc_samples.RDS')
sm <- readRDS('/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/snpStats/t1d_wtccc_bootstrap_USE.RDS')

## for time being just sample w/o stratifiying by b58c region

getSample <- function(DT,status,n){
    idx=which(DT$t1d==status)
    regions=split(idx,DT$b58region[idx])
    sidx <- do.call('c',lapply(regions,function(x){
      nr <- ceiling((length(x)/length(idx) *n))
      sample(x,nr)
    }))
    ## likely more than n so downsample accordingly
    sort(sample(sidx,n))
}

## next fit log regressions and get summary stats for projection
sample.DT<-data.table(samples)
gw <- lapply(c(500,1000,1900),function(nm){
  message(sprintf("Processing %d",nm))
  idx <- c(getSample(sample.DT,1,nm),getSample(sample.DT,2,nm))
  bs.samples <- samples[idx,]
  bs.snps <- sm[idx,]
  n1 <- sum(bs.samples$cc==1)
  n <- nrow(bs.samples)
  rownames(bs.snps) <- rownames(bs.samples)
  res <- snp.rhs.estimates(cc ~ b58region,family='binomial',data=bs.samples,snp.data=bs.snps,uncertain = TRUE)
  or <- sapply(res,'[[','beta')
  varb <- sapply(res,'[[','Var.beta')
  p <- sapply(res,function(x) 2*pnorm(abs(x$beta/sqrt(x$Var.beta)),lower.tail=FALSE))
  z<-sapply(res,function(x) x$beta/sqrt(x$Var.beta))
  data.table(pid=names(res),or=exp(or),Z=z,p.value=p,trait=sprintf("T1D_%d_%d",nm,nm),se.b=sqrt(varb),n=nm,n1=nm)
})

saveRDS(gw,file="/home/ob219/rds/hpc-work/as_basis/callibration/wtccc_t1d_gwas_samples/all_wtccc.RDS")



## next we load in basis file and do the projection

BASIS.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/basis_ws.RDS"
SHRINK.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/shrink_ws.RDS"
SUPPORT.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/snp_support_bootstrap_USE.tab"
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_bootstrap/'
basis <- readRDS(BASIS.FILE)
shrink <- readRDS(SHRINK.FILE)
support <- fread(SUPPORT.FILE)
setkey(support,pid)

w.DT <- data.table(pid=rownames(basis$rotation),basis$rotation)

var <- lapply(gw,function(gw.DT){
  message(sprintf("Processing %d",unique(gw.DT$n1)))
  gw.DT[,emp_se:=se.b]
  gw.DT[,c('chr','position'):=tstrsplit(pid,':')]
  setkey(gw.DT,pid)
  gw.DT <- gw.DT[support[,.(pid,ld.block)]]
  compute_proj_var(gw.DT,w.DT,shrink,DEFAULT.SNPSTATS.DIR,'ws_emp',quiet=TRUE)
})

var.DT<-data.table(cases=c(500,1000,1900),total=2*c(500,1000,1900),do.call('rbind',var))
var.DT <- melt(var.DT,id.vars=c('cases','total'),measured.vars=paste('PC',1:10,sep=""))

saveRDS(var.DT,file="~/tmp/act_var_bootstrap.RDS")
