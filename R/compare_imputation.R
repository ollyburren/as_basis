library(data.table)
data.dir<-'/Users/oliver/DATA/AS_BASIS/t1d.imp.comparison';
(load(file.path(data.dir,'t1d_ncooper.RData')))
(load(file.path(data.dir,'t1d_wallace.RData')))
cooper<-ndat[,c('id','chromosome','position','a0','a1','beta.meta','se.meta'),with=FALSE]


# compare the OR across all chromosomes
dat<-rbindlist(dat)
wallace<-dat[,c('rsid','chr','position','alleleA','alleleB','beta.meta','se.meta'),with=FALSE]
setnames(wallace,names(cooper))

to.change<-c('id','a0','a1','beta.meta','se.meta')
setnames(cooper,to.change,paste('cooper',to.change,sep='.'))
setnames(wallace,to.change,paste('wallace',to.change,sep='.'))
wallace$chromosome<-as.numeric(wallace$chromosome)
setkeyv(wallace,c('chromosome','position'))
setkeyv(cooper,c('chromosome','position'))
## filter where we have no positional information as cannot do anything with this
wallace<-subset(wallace,!is.na(chromosome))
cooper<-subset(cooper,!is.na(chromosome))

##merge

m<-wallace[cooper]
## there are many instances where wallace not found in cooper (possible MAF filtering)
length(idx<-which(is.na(m$wallace.a0)))
## if we look at where they do overlap
m<-m[-idx,]
## remove where the alleles don't agree
length(idx<-which(m$wallace.a0!=m$cooper.a0))
m<-m[-idx,]

## compute pvals

m$cooper.p<-2*pnorm(abs(m$cooper.beta.meta/m$cooper.se.meta),lower.tail=FALSE)
m$wallace.p<-2*pnorm(abs(m$wallace.beta.meta/m$wallace.se.meta),lower.tail=FALSE)

## compare beta

## do the signs always agree ?

length(idx<-which(sign(m$wallace.beta.meta) != sign(m$cooper.beta.meta)))

## how many are GWAS - sig 

diff.dir<-m[idx,]
odd<-subset(diff.dir,(cooper.p<5e-8 | wallace.p<5e-8) & nchar(wallace.a0)==1 & nchar(wallace.a1)==1 & nchar(cooper.a0)==1 & nchar(cooper.a1)==1 )

library(ggplot2)

ggplot(m,aes(x=wallace.beta.meta,y=cooper.beta.meta)) + geom_hex()

## what is going on with beta between the two studies.

