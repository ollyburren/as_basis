## diagram to explain approximate BayesFactors
library(Gviz)
library(data.table)
library(GenomicRanges)
library(biomaRt)

## get the data from /scratch/ob219/pid/GWAS

ra<-readRDS("RA.RDS")
t1d<-readRDS("T1D.RDS")

## attempt to show details for CLEC16A region

## find out LD blocks

#chr2:204732511 204738683
#chr21:43823971-43867790

#UBASH3A
rchr<-'21'
rstart<-43823971
rend<-43867790

#GLIS3 3898646-3901248
rchr <- '9'
rstart<-3824128
rend<-4300035

#IL2
rchr <- '4'
rstart<-123372626
rend<-123377650

#CTLA4
rchr<-'2'
rstart<-204732511
rend<-204738683

#PADI4
rchr<-'1'
rstart<-17634690
rend<-17690495

ld.block<-unique(subset(t1d,chr==rchr & between(position,rstart,rend))$ld)
ld.blocks <- c(ld.block-1,ld.block,ld.block+1)
ld.blocks<-ld.block

ra.dat <- subset(ra,ld %in% ld.blocks & !duplicated(pid))
t1d.dat <- subset(t1d,ld %in% ld.blocks & !duplicated(pid))

pid.in <-intersect(ra.dat$pid,t1d.dat$pid)

ra.dat <- ra.dat[pid %in% pid.in ,]
t1d.dat <- t1d.dat[pid %in% pid.in ,]

ra.dat[,trait:='RA']
t1d.dat[,trait:='T1D']


comb.dat<-rbind(ra.dat,t1d.dat)
comb.dat[,mlp:=-log10(p.val)]
comb.dat[,lp0:=log(1-ppi)]

## take the q_i from genomewide calcs at 0.001
emp<-0.001
tmp<-comb.dat[,list(q_i=1-exp(sum(lp0,na.rm=TRUE))),by=c('ld','pid')]
po<-emp/(1-emp)
## note that emp is for h1 that beta != 0 therefore we need to take reciprocal as equation assumes pi_0 - see notes
po<-1/po
tmp[,uABF:=po*(q_i/(1-q_i))]
## set an upper limit to BF (here its upper 0.0001 percentile)
BIG<-quantile(tmp[is.finite(tmp$uABF),]$uABF,prob=0.9999)
tmp[is.infinite(uABF) | uABF > BIG, uABF:=BIG]


logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

basis_pp<-function(BF,emp_pi){
  lABF<-log(BF)
  tABF <- c(lABF,0)
  vpi_i<-c(rep(emp_pi,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(emp_pi)-sBF)
}

tmp<-tmp[,bshrink:=basis_pp(uABF,emp),by=ld][,.(bshrink),by=pid]
setkey(tmp,pid)



m.dat<-melt(comb.dat,id.vars=c('pid','chr','position','ld','trait'),measure.vars=c('ppi','mlp'))
comb.dat <- dcast(m.dat,chr+position+ld+pid~trait+variable)

setkey(comb.dat,pid)

comb.dat<-comb.dat[tmp]



comb.gr <- with(comb.dat,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L),RA.ppi=RA_ppi,T1D.ppi=T1D_ppi,RA.mlp=RA_mlp,T1D.mlp=T1D_mlp,shrinkage=bshrink))


e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")



## attempt to combine



ideo<-fread("/home/ob219/scratch/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

tracks<-list()

tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(comb.gr),bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        chromosome = seqlevels(comb.gr), start = min(t1d.dat$position), end = max(t1d.dat$position),
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T))
displayPars(tracks$gene)$background.title <- 'black'
tracks$combmlp<- DataTrack(name="-log10(P)",comb.gr[,c('RA.mlp','T1D.mlp')],groups=c('RA','T1D'))
displayPars(tracks$combmlp)$background.title <- 'black'
tracks$rapost<- DataTrack(name="Posterior RA",comb.gr,data=comb.gr$RA.ppi)
displayPars(tracks$rapost)$background.title <- 'dodgerblue'
tracks$t1dpost<- DataTrack(name="Posterior T1D",comb.gr,data=comb.gr$T1D.ppi,col="magenta")
displayPars(tracks$t1dpost)$background.title <- 'magenta'
tracks$shrinkage<- DataTrack(name="Combined Weight",comb.gr,data=comb.gr$shrinkage,col="purple")
displayPars(tracks$shrinkage)$background.title <- 'purple'



plotTracks(tracks,from=min(t1d.dat$position),to=max(t1d.dat$position),sizes=c(1,1,1,1,1,1,3))
