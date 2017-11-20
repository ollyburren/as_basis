## ok given a genotype compute the expected posterior log(OR)
#library(devtools)
#install_github('ollyburren/cupcake')

library(optparse)

option_list = list(
        make_option(c("-o", "--out"), type="character", default=NULL,
              help="output dir", metavar="character"),
        make_option(c("-c", "--chromosome"), type="character",
                    help="chromosome to parse", metavar="character"),
        make_option(c("-s", "--samplefile"), type="character", default=NULL,
                                help="List of samples to process", metavar="character")
        )


## This script processes RAW summary stats it uses a pregenerated list of
## here we define the posterior lookup - this needs to be calculated for a given basis
## see lor_lookup.R for more details
support.dir<-'/scratch/ob219/as_basis/support_tab'
lor.lu.file <- file.path(support.dir,'lor_posterior.tab')
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
#message(args)

library(cupcake)

## read in samples file
if(FALSE){
args <- list(
  samplefile = '/home/ob219/scratch/as_basis/jia_ind_analysis/splits/xaa',
  chr = '1',
  out = '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_basis5'
)
}
samples <- scan(file=args$samplefile,"character()")

dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
setnames(dat,c('f','lor.00','lor.01','lor.11'))
## lets create a model for each of these
colnames(dat) <- make.unique(colnames(dat))
m <- reshape2::melt(dat,"f")


#library(ggplot2)
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
ref_af_file<-file.path(support.dir,'as_basis_snps_with_af.tab')
#basis.snps <- fread(ref_af_file)[,pid:=paste(chr,position,sep=':')]
basis.snps <- fread(ref_af_file)
setkey(basis.snps,pid)

# really we need to know which allele goes with which otherwise this won't work
# need a different support file

#(load("all_EUR_support.RData"))
# subset to just as.basis snps
#all.eur[,pid:=paste(chr,position,sep=':')]
#setkey(all.eur,pid)
#all.eur.basis <- subset(all.eur,pid %in% basis.snps$pid)
#write.table(all.eur.basis,file='/scratch/ob219/as_basis/support_tab/as_basis_snps_with_af.tab',quote=FALSE,row.names=FALSE,sep="\t")

DTl <- split(basis.snps,basis.snps$chr)
DTl <-  lapply(DTl,setkey,key='pid')
library(annotSnpStats)


## loop over each chromosome and create a matrix of expected log OR
 target.chr <- args$chr
 gt.file <- file.path('/scratch/ob219/as_basis/JIA_basis_annotSnpStats',sprintf('annotsnpstats-%s.RData',target.chr))
 message(sprintf("Loading %s",gt.file))
 X <- get(load(gt.file))

 gwas.snps <- data.table(snps(X))[,c('pid','order'):=list(paste(chromosome,position,sep=':'),1:.N),]
 ## only care about controls
 sample.idx <- which(samples(X)$ID_1 %in% samples)
 if(length(sample.idx) != length(samples))
  stop("Missing some some samples in gt file !")
 message("Converting to snpMatrix object")
 sm<-as(X,"SnpMatrix")

 sm <- sm[sample.idx,]
 ## This is super slow perhaps convert to a snpMatrix object ?
 #sm <- as(X[sample.idx,snp.idx],"SnpMatrix")
 #For each SNP get lor for each snp config
 #gwas.snps <- gwas.snps[order %in% snp.idx,]

 gwas.snps <- gwas.snps[,AF:=af.wrt.a2]
 gwas.snps[,AF:=round(af.wrt.a2,digits=3)]
 ## avoid asymptotics
 gwas.snps[AF<0.01,AF:=0.01]
 gwas.snps[AF>0.99,AF:=0.99]
 #note we use AF.a2
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
 if(any(is.na(Xst$lor))){
   stop("Things don't join up aborting")
 }
 #Xst[is.na(lor),]
 Xst <- Xst[,.(sample,pid,lor)]

 ref.DT <- subset(basis.snps,chr==target.chr)
 DT<-gwas.snps[,.(chromosome,position,allele.1,allele.2)]
 setnames(DT,c('chr','position','a1','a2'))

 flip.idx <- flip_allele(DT,ref.DT)
 Xst[,c('or','flipped'):=list(exp(lor),0)]
 Xst[pid %in% DT[flip.idx,]$pid,c('or','flipped'):=list(exp(lor*-1),1)]
 of <- sprintf("chr%s_%s.RDS",target.chr,basename(args$samplefile))
 message(sprintf("Writing results to %s",of))
 saveRDS(Xst[,.(sample,pid,or,flipped)],file=file.path(args$out,of))








 # also store the support file
