## ok given a genotype compute the expected posterior log(OR)

library(cupcake)
support.dir<-'/scratch/ob219/as_basis/support_tab'
lor.lu.file <- file.path(support.dir,'lor_posterior.tab')
dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
setnames(dat,c('f','lor.00','lor.01','lor.11'))
## lets create a model for each of these

library(ggplot2)
colnames(dat) <- make.unique(colnames(dat))
m <- reshape2::melt(dat,"f")
#ggplot(m[grep("%", m$variable,invert=TRUE),],
#       aes(x=f,y=value,col=variable,group=variable)) + geom_point() + geom_smooth() + labs(x="AF",y="log OR")


gtn<-c('00','01','11')

mods<-lapply(paste('lor',gtn,sep='.'),function(x){
  form <- sprintf('%s ~ f',x)
  loess(form,dat)
})

names(mods) <- gtn

## next load in some genotype data

gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
ref_af_file<-file.path(support.dir,'as_basis_snps.tab')
basis.snps <- fread(ref_af_file)[,pid:=paste(chr,position,sep=':')]
DTl <- split(basis.snps,basis.snps$chr)
library(annotSnpStats)

## loop over each chromosome and create a matrix of expected log OR
  chr <- 1
  gt.file <- file.path('/scratch/wallace/JIA-2017-data/',sprintf('annotsnpstats-%d.RData',chr))
  X <- get(load(gt.file))

 gwas.snps <- data.table(snps(X))[,c('pid','order'):=list(paste(chromosome,position,sep=':'),1:.N),]
 snp.idx <- which(gwas.snps$pid %in% DTl[[chr]]$pid)
 if(length(snp.idx) != nrow(DTl[[chr]]))
  stop("Not all SNPs in the basis")
 ## only care about controls
 sample.idx <- which(samples(X)$phenotype==1)
 sm<-as(X,"SnpMatrix")
 sm <- sm[sample.idx,snp.idx]
 ## This is super slow perhaps convert to a snpMatrix object ?
 #sm <- as(X[sample.idx,snp.idx],"SnpMatrix")
 #For each SNP get lor for each snp config
 gwas.snps <- gwas.snps[order %in% snp.idx,]
 gwas.snps <- gwas.snps[,c('lor.00','lor.01','lor.11'):=lapply(mods,predict,AF)]
 Xs <- data.table((apply(matrix(sm,nrow(sm),ncol=ncol(sm)),2,as.character)))
 setnames(Xs,colnames(sm))
 Xs[,sample:=rownames(sm)]
 Xst <- melt(Xs,id.vars='sample')
 setkey(Xst,variable)
 setkey(gwas.snps,ID)
 Xst <- Xst[gwas.snps[,.(ID,pid,lor.00,lor.01,lor.11)]]
 Xst[value=='01',lor:=lor.00]
 Xst[value=='02',lor:=lor.01]
 Xst[value=='03',lor:=lor.11]
 Xst[is.na(lor),]
 Xst <- Xst[,.(sample,pid,lor)]

 # also store the support file
