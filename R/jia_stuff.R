library(data.table)

# read in the data

setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")

DT<-getGWASData()

## read in the basis as we can only accommodate SNPs that are in here so can filter out the rest

cc.jia<-readRDS("/home/ob219/scratch/jia/by.trait/cc.RDS")
cc.jia[,id:=paste(chr,position,sep=':')]
setkey(cc.jia,id)
filt<-cc.jia[cc.jia$id %in% unique(DT$id),]
## construct DT that is compatible with basis
skel<-filt[,.(chr,position,id,a1,a2,beta,p.val)]
## so beta are wrt to a2 but we need to line these up with DT
alleles<-unique(DT[,.(id,a1,a2,name,risk.allele.freq)])
setnames(skel,c('a1','a2'),c('sa1','sa2'))
skel<-skel[alleles]

## match

midx<-with(skel,which(sa1==a1 & sa2==a2))
sidx<-with(skel,which(sa1==a2 & sa2==a1))
nomidx<-setdiff(1:nrow(skel),c(midx,sidx))
