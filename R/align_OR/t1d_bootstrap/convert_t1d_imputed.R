library(snpStats)
library(optparse)
library(data.table)


TEST <- FALSE
MAF.THRESH <- 0.0001
CALL.RATE <- 0.1
HWE.Z.THRESH <- 20

option_list = list(
        make_option(c("-i", "--in_file"), type="character", default=NULL,
                help="GWAS file to filter", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Directory to write results", metavar="character"),
        make_option(c("-s", "--sample_file_dir"), type="character", default=NULL,
                help="sample file dir", metavar="character")
        )


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if(!file.exists(args$in_file))
    stop("Can't find source data files")
}else{
  #for testing
  message("IN TESTING MODE ======>!")
  args <- list(
    in_file='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/imputed-wtccc-22.gen.gz',
    sample_file_dir='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/',
    out_dir='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/snpStats'
  )
}

study <- gsub("imputed\\-([^\\-]+)\\-([^\\.]+)\\.gen\\.gz","\\1",basename(args$in_file))
chr <- gsub(".*-([^\\.]+)\\.gen\\.gz","\\1",basename(args$in_file))

## use this to construct sample file
sample.file <- file.path(args$sample_file_dir,sprintf("samples-wtccc-%s.gz",chr))
if(!file.exists(sample.file))
  stop(sprintf("Sample file %s does not exist",sample.file))
samples <- fread(paste("zcat",sample.file),skip=1L)
setnames(samples,c('sampleid','t1d','sex','collection','caucasian','cc','id','pedigree','mother','father','cohort'))
s <- read.impute(args$in_file)
## remove duplicated SNP names
dup.idx <- which(duplicated(colnames(s)))
if(length(dup.idx)>0){
  s<-s[,-dup.idx]
}

## should double check with input into impute.
imp.sample.file <- file.path(dirname(args$in_file),sprintf("imputed-%s-%s.gen_filt.sample",study,chr))
if(file.exists(imp.sample.file)){
  s2<-fread(imp.sample.file,sep=" ")[ID_1!=0,]
  if(sum(s2$ID_1==samples$sampleid)!=length(samples$sampleid)){
    message("Sample mismatch between imputed and input data")
    stop()
  }
}else{
  message(sprintf("No impute sample %s cannot check",imp.sample.file))
}

## now we compute MAF and remove those that are below threshold or are not in HWE
idx<-which(samples$t1d==1)
sum <- col.summary(s[idx,])
snp.keep <- with(sum,which(Call.rate>CALL.RATE & MAF>MAF.THRESH & abs(z.HWE)<HWE.Z.THRESH))
s.filt <- s[,snp.keep]
fname <- gsub(".gen.gz",".RData",basename(args$in_file))
snps <-data.table(id=colnames(s.filt))
snps[,c('rs','position','a1','a2'):=tstrsplit(id,':')]
snps[,chr:=chr]
obj <- list(sm=s.filt,snps=snps,samples=samples)
save(obj,file=file.path(args$out_dir,fname))
message(paste("Success",file.path(args$out_dir,fname)))

if(FALSE){
  ## get a  list of SNPs for wtccc portion
  library(data.table)
  library(snpStats)
  files <- files <- list.files(path='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/tmp/snpStats',pattern="*.RData",full.names=TRUE)
  all.obj <- lapply(files,function(x) get(load(x)))
  sm<-do.call('cbind',lapply(all.obj,'[[','sm'))
  snps<-rbindlist(lapply(all.obj,'[[','snps'))
  samples<-all.obj[[1]]$samples
  rownames(sm) <- samples$sampleid
  ## add in UK regions
  ofs <- data.table(get((load("/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-gwas-data/corrected-support/WTCCC-sample-support.RData"))))
  lu<-split(ofs$b58region,ofs$affy.id)
  samples[,b58region:=unlist(lu[sampleid])]
  ## generate indices for bootstrapping
  library(magrittr)
  library(snpStats)

  ## chris' bootstrap function
  getBoot <- function(sample.DT,status,n.boot){
    idx=which(sample.DT$t1d==status)
    regions=split(idx,sample.DT$b58region[idx])
    replicate(n.boot, sapply(regions, function(x){
      if(length(x)==1){
        return(x)
      }else{
        sample(x,replace=TRUE)
      }
    }) %>% unlist(.,use.names=FALSE))
  }

  n.cases<-100
  n.controls<-100
  n.boostraps<-500
  ## create bootstrap objects as quicker to do this first and then index
  idx=c(sample(which(samples$t1d==2 & samples$b58region!=-1),n.cases),
  sample(which(samples$t1d==1 & samples$b58region!=-1),n.controls))
  bs.sm <- sm[idx,]
  bs.samples <- samples[idx,]
  bs.samples[,cc:=t1d-1]
  #each column is a bootstrap sample
  bs.mat <- rbind(getBoot(bs.samples,1,n.boostraps),getBoot(bs.samples,2,n.boostraps))
  ## fit logistic
  bs.samples <- as.data.frame(bs.samples)
  rownames(bs.samples) <- bs.samples$sampleid


  ##test code for how to do bootstrap GWAS - run on queue eventually does for 20 SNPs
  test.boot <- bs.mat[,1]
  sample.dat <- bs.samples[test.boot,]
  snp.dat <- as(bs.sm[test.boot,],"SnpMatrix")
  rownames(snp.dat) <- rownames(sample.dat)
  test=snp.rhs.estimates(cc ~ b58region,family='binomial',data=sample.dat,snp.data=snp.dat,sets=1:20)

}
