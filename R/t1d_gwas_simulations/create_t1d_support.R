library(snpStats)
library(parallel)
library(magrittr)
library(rmeta)
library(annotSnpStats)

d <- "/scratch/wallace/t1d-barrett-gwas-data/corrected-support"
(load(file.path(d,"T1DGC-sample-support.RData")))
(load(file.path(d,"WTCCC-sample-support.RData")))

library(snpStats)
dg <- file.path(d,"../snpStats-genotypes")
files <- list.files(dg,pattern="WTCCC",full=TRUE)

## output
dr <- '/home/ob219/scratch/as_basis/t1d_gwas_sim/support/'

for(f in files) {
message(f)
chr <- gsub(".*WTCCC.|.clean.*","",f)
message(sprintf("Processing %s",chr))
outf <- file.path(dr,paste0("results-",chr,".tsv"))
if(file.exists(outf)){
  next
}
wtccc <- eval(as.symbol(load(f)[1]))
t1dgc <- eval(as.symbol(load(sub("WTCCC","T1DGC",f))[1]))
c58 <- eval(as.symbol(load(sub("WTCCC","C58",f))))
sample.dups <- subset(wtccc.sample.support,affy.id %in% rownames(c58) & sanger.id %in% rownames(t1dgc))
ss <- eval(as.symbol(load(file.path(d,paste0("snp-support-",chr,".RData")))))
cols <- intersect(colnames(wtccc.sample.support),colnames(t1dgc.sample.support))
m1 <- match(rownames(wtccc),wtccc.sample.support$affy.id)
wtccc.support <- wtccc.sample.support[m1[!is.na(m1)],cols]
wtccc.snps <- wtccc[which(!is.na(m1)),]
m1 <- match(rownames(wtccc),wtccc.sample.support$sanger.id)
m2 <- match(rownames(t1dgc),t1dgc.sample.support$sample)
t1dgc.support <- rbind(wtccc.sample.support[m1[!is.na(m1)],cols],t1dgc.sample.support[m2[!is.na(m2)],cols])
t1dgc.snps <- t1dgc[c(which(!is.na(m1)),which(!is.na(m2))),]
## save lined up data
wtccc.support$pedigree <- wtccc.support$id <- rownames(wtccc.snps)
wtccc.support$father <- wtccc.support$mother <- 0
t1dgc.support$pedigree <- t1dgc.support$id <- rownames(t1dgc.snps)
t1dgc.support$father <- t1dgc.support$mother <- 0
rownames(wtccc.support) <- rownames(wtccc.snps)
rownames(t1dgc.support) <- rownames(t1dgc.snps)
ss$wtccc.1 <- substr(ss$wtccc.affy_alleles,1,1)
ss$wtccc.2 <- substr(ss$wtccc.affy_alleles,3,3)
ss$t1dgc.1 <- substr(ss$ilmn_alleles,1,1)

## switch alleles

mat <- matrix(FALSE,nrow(ss),5,dimnames=list(NULL,c("NA","nochange","rev","comp","revcomp")))
x <- ss$wtccc.affy_alleles
y <- ss$ilmn_alleles
mat[,"NA"] <- x=="" | y==""
use <- !mat[,"NA"]
## nochange
mat[use , "nochange" ] <- x[use]==y[use]
mat[use, "rev"] <- x[use]==g.rev(y[use])
mat[use,"comp"] <- x[use]==g.complement(y[use])
mat[use,"revcomp"] <- x[use]==g.rev(g.complement(y[use]))
rowMeans(mat) %>% table() ## should all be 0.2 - ie only 1 column selected
sw <- mat[,"rev"] | mat[,"revcomp"]
## switch illumina
t1dgc.snps <- switch.alleles(t1dgc.snps,colnames(t1dgc.snps) %in% rownames(ss)[sw])
ss$i.switched <- FALSE
ss$i.switched[sw] <- TRUE

## basic summary
cs.a0 <- col.summary(wtccc.snps[wtccc.support$t1d==0,])
cs.a1 <- col.summary(wtccc.snps[wtccc.support$t1d==1,])
cs.i0 <- col.summary(t1dgc.snps[t1dgc.support$t1d==0,])
cs.i1 <- col.summary(t1dgc.snps[t1dgc.support$t1d==1,])
cs.a <- merge(cs.a0[,c("RAF"),drop=FALSE],cs.a1[,c("RAF"),drop=FALSE],
              by=0,suffixes=c(".ctl",".cse"))
cs.i <- merge(cs.i0[,c("RAF"),drop=FALSE],cs.i1[,c("RAF"),drop=FALSE],
              by=0,suffixes=c(".ctl",".cse"))
cs.a$Row.names %<>% as.character()
cs.i$Row.names %<>% as.character()
cs <- merge(cs.a,cs.i,by="Row.names",all=TRUE,suffixes=c(".aff",".ill"))
pdf(file.path(dr,paste0(chr,'QC.pdf')))
plot(cs$RAF.ctl.aff,cs$RAF.ctl.ill)
dev.off() # check

  ## make imputation rules
AFFY <- c58[as.character(sample.dups$affy.id),]
ILL <- t1dgc[as.character(sample.dups$sanger.id),]
A.unique <- setdiff(colnames(AFFY),colnames(ILL))
I.unique <- setdiff(colnames(ILL),colnames(AFFY))

a2i <- snp.imputation(X=AFFY, Y=ILL[,I.unique], pos.X=ss[colnames(AFFY),"b36_start"], pos.Y=ss[I.unique,"b36_start"])
i2a <- snp.imputation(Y=AFFY[,A.unique], X=ILL, pos.Y=ss[A.unique,"b36_start"], pos.X=ss[colnames(ILL),"b36_start"])
save(list=c('wtccc.support','t1dgc.support','wtccc.snps','t1dgc.snps','ss','a2i','i2a'),
  file=file.path(dr,paste0(chr,"_support.RData")))
}
