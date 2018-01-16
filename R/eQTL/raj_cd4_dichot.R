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

#load in filtered genotype data
raj_cd4 <- readRDS("/scratch/ob219/as_basis/raj_cd4_ind_analysis//raj-cd4-basis_genotypes.RDS")
# load in the expression data
translate <- readRDS('/scratch/wallace/twas/model_output/trans_raj-cd4.rds')


load("/scratch/wallace/twas/raj-cd4-expression.RData")
expr.DT<-data.table(t(expr))

## there are some duplicates in translate file - checked these out and appears that happens on Raj side
## so we can ignore.

translate <- translate[names(translate) %in% colnames(expr)]

lu.s <- split(names(translate),translate)
s.DT<-data.table(samples(raj_cd4))
s.DT[,chip:=unlist(lu.s[pedigree])]


## load in the annotation file as no point in checking things that have duplicate probes

raj.anno <- data.table(readRDS("/scratch/wallace/twas/gene_details.rds"))
duplicated_probes<-raj.anno[duplicated(raj.anno$affy_hugene_1_0_st_v1),]$affy_hugene_1_0_st_v1
fexpr.DT <- expr.DT[,-which(names(expr.DT) %in% duplicated_probes),with=FALSE]

## next we want to create a set of phenotype files for whether expression value is above or below mean expression
all.mean.expression<-apply(fexpr.DT,2,mean)
## in practice the mean seems to be zero so some standardisation already in place ?

outcome.mean <- ((t(fexpr.DT)  > colMeans(fexpr.DT)) %>% t ) * 1
outcome.median <- apply(fexpr.DT,2,function(z) z > median(z)) * 1

rownames(outcome.mean) <- colnames(expr)
pheno.mean <- cbind(s.DT,outcome.mean[s.DT$chip,])

rownames(outcome.median) <- colnames(expr)
pheno.median <- cbind(s.DT,outcome.median[s.DT$chip,])

## computes OR w/o fitting logistic function so very quick
## Y is a matrix or vector of outcomes {0,1}
## G is a matix or vector  of genotypes {0,1,2}
computeORMat <- function(Y,G,n0=214,n1=212){
  # ploidy is two - two groups are not identical if we use median
  if(is.matrix(Y) & is.matrix(G))
    stop("Both X and Y are matrices in computeORMat")
  unx <- colSums((1-Y) * G)
  ex <- colSums(Y * G)
  #(log(n0-unx) + log(ex)) - (log(unx) + log(n1-ex)) ##==
  log(n0-unx) + log(ex) - log(unx) - log(n1-ex)
}

computeSEORMat <- function(Y,G,n0=214,n1=212){
  if(is.matrix(Y) & is.matrix(G))
    stop("Both X and Y are matrices in computeORMat")
  unx <- colSums((1-Y) * G)
  ex <- colSums(Y * G)
  a <- (n0-unx)
  b <- unx
  c <- (n1-ex)
  d <- ex
  sqrt(1/(n0+n1)) * sqrt(1/a + 1/b + 1/c + 1/d)
}

bsnps <- fread("/scratch/ob219/as_basis/support_tab/as_basis_snp_support_eQTL.tab")
anno <- data.table(snps(raj_cd4))[,rorder:=1:.N]
setkey(anno,pid)
setkey(bsnps,pid)
sdt <- anno[bsnps][order(rorder),]
setnames(sdt,c('allele.1','allele.2'),c('a1','a2'))
sdt[,flag:='unprocessed']
sdt[a1==ref_a1 & a2==ref_a2,flag:='match']
sdt[a1==ref_a2 & a2==ref_a1,flag:='flip']
sdt[a1==comp(ref_a1) & a2==comp(ref_a2),flag:='match_rc']
sdt[a1==comp(ref_a2) & a2==comp(ref_a1),flag:='flip_rc']

## align alleles with the basis
sm <- as(raj_cd4,'SnpMatrix')
snp.ids <- sdt$pid
saveRDS(snp.ids,file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/snps.RDS")
#sm<-matrix(sm,nrow=nrow(sm),ncol=ncol(sm),dimnames=list(NULL,sdt$pid))
sm<-matrix(sm,nrow=nrow(sm),ncol=ncol(sm))
fix <- which(sdt$flag %in% c('flip','flip_rc'))
g.mat <- apply(sm,2,as.numeric)
g.mat[,fix] <- abs(3-g.mat[,fix])
g.mat[,-fix] <- g.mat[,-fix] -1
saveRDS(g.mat,file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/gt_matrix.RDS")

ot.mat <- as.matrix(pheno.median[,8:ncol(pheno.median)])
saveRDS(ot.mat,file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/outcome_matrix.RDS")

# like this to do each SNP vs one probe
all.snps.vs.one.probe <- computeORMat(ot.mat[,1],g.mat)
# like this to each probe vs a SNP
all.probes.vs.one.snp <- computeORMat(ot.mat,g.mat[,1])

## practice with 10 probes

#system.time(res<-lapply(1:10,function(pidx){
#  pout <- ot.mat[,pidx]
#  computeORMat(pout,g.mat)
#}))

## we then need to package this and project onto basis !

library(cupcake)
support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support_eQTL.tab')
m_file<-file.path(support.dir,'as_manifest_december.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'

basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/shrink.RDS")
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
## save this
saveRDS(pc.emp,file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/pc.RDS")


## this works but is too slow - really we need to work out how to do on the Q
library(parallel)
all <- mclapply(1:1000,function(j){
  message(j)
  # flip alleles and adjust genotype dosage so {0,1,2} not {1,2,3}
  if(sdt$flag[j] %in% c('flip','flip_rc')){
    g <- abs(3-as.numeric(sm[,j]))
  }else{
    g <- as.numeric(sm[,j]) -1
  }
  computeORMat(ot.mat,g)
  #sapply(8:ncol(pheno.median),function(p) computeOR(pheno.median[[p]],g,length(g)))
},mc.cores=1)


## code an plots to see that the computation above works
# there is close but not exact agreement between logistic and 2 x 2 fit
# I wonder if this is a function of the minor allele freq as this will
# affect the standard errors and hence the logistic fit.

if(FALSE){
  sm <- as(raj_cd4,'SnpMatrix')
  sm<-matrix(sm,nrow=nrow(sm),ncol=ncol(sm))
  j<-1
  g <- abs(as.numeric(sm[,j])) - 1
  ## check OR calc for a few probes
  i<-1:10000 + 7
  lors <- lapply(i,function(i){
    clor <- computeOR(pheno.median[[i]],g,length(g))
    dat <- data.table(Y=pheno.median[[i]],gt=g)
    setnames(dat,c('Y','gt'))
    glor <- coefficients(summary(glm(Y~gt,data=dat,family=binomial())))[2,1]
    data.table(clor,glor)
  })
  library(ggplot2)
  ggplot(rbindlist(lors),aes(x=abs(clor),y=abs(glor))) + geom_point() + geom_abline(intercept=0,slope=1,color='red') +
  xlab("log(OR) from 2 x 2") + ylab("log(OR) from logistic")
}

library(magrittr)
ret <- lapply(1:10,function(pidx){
  pout <- ot.mat[,pidx]
  computeORMat(pout,g.mat)
}) %>% do.call("rbind",.)

## as input expect OR

ret <- exp(ret) %>% t %>% data.table
ret[,snp:=snp.ids]
setnames(ret,c(as.character(colnames(ot.mat)[1:10]),'pid'))
DT <- melt(ret,id.vars='pid')
setnames(DT,c('pid','sample','or'))
## some OR are inf assume because missing allele ?
DT[is.infinite(DT$or) | DT$or==0 ,or:=1]

setcolorder(DT,c('sample','pid','or'))
setkey(DT,pid)
setnames(DT,'sample','trait')
mat.emp <- create_ds_matrix(DT,shrink.DT,'emp')
proj<-predict(pc.emp,newdata=mat.emp)

## to check stuff we want a random set of OR for this we need the standard error of OR

## first we compute this for a random probeset - doesn't matter which here
## we just take the first one
all.snps.vs.one.probe.SE <- computeSEORMat(ot.mat[,1],g.mat)
