## This script simulates what happens under the null i.e when the scaled beta's i.e. gamma hats are equal to one another.

## assumes that cupcake has been loaded with (devtools)
library(cupcake)
library(optparse)

TEST <- FALSE # set to true when debugging
DEFAULT.SUPPORT.DIR <- '/home/ob219/scratch/as_basis/support_tab'
DEFAULT.SNPSTATS.DIR <- '/home/ob219/scratch/as_basis/snpStats/basis_1kg/'
DEFAULT.N.SIMS <- 50
DEFAULT.GWAS.DATA.DIR <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
#DEFAULT.MANIFEST.FILE <- 'as_basis_manifest_with_jia_cc.tab'
DEFAULT.MANIFEST.FILE <- 'as_manifest_february_2018.tab'
DEFAULT.LD.BLOCK.FILE <- 'all.1cM.tab'
DEFAULT.REF.AF.FILE <- 'as_basis_snp_support.tab'


## this routine takes a set of vcf files from a reference panel say 1KG and extracts the variants that are used in our basis.
#In this form it won't work (needs basis info for filtering but acts as a reference)

filterVCF <- function(){
  bcft <- '~/bin/bcftools-1.4/bcftools'
  vcf.dir<-'/home/ob219/scratch/DATA/1kgenome/VCF/EUR/by.chr.phase1/'
  out.dir <- '/home/ob219/scratch/as_basis/snpStats/basis_1kg/'
  s.DT <- split(study1.DT[,.(chr,position)],study1.DT$chr)
  for(chr in names(s.DT)){
      message(sprintf("Processing %s",chr))
      tmp<-s.DT[[chr]][order(as.numeric(chr),as.numeric(position)),]
      tmp.file <- tempfile(pattern = "vcf_tmp", tmpdir = tempdir())
      write.table(tmp,file=tmp.file,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      #vcf.file <- file.path(vcf.dir,sprintf("ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz",chr))
      vcf.file <- file.path(vcf.dir,sprintf("chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz",chr))
      obj <- vcf2snpmatrix(vcf.file,bcft,tmp.file,TRUE)
      unlink(tmp.file)
      out.file <- file.path(out.dir,sprintf("%s_1kg.RData",chr))
      save(obj,file=out.file)
    }
}

## to prevent unneccsary computation we can compute the basis ahead of time and just read that in

cacheBasis <- function(args){
  support.dir <- args$support_dir
  ref_af_file<-file.path(support.dir,DEFAULT.REF.AF.FILE)
  ld_file<-file.path(support.dir,DEFAULT.LD.BLOCK.FILE)
  m_file<-file.path(support.dir,DEFAULT.MANIFEST.FILE)
  basis.DT<-get_gwas_data(m_file,ref_af_file,args$gwas_data_dir)
  shrink.DT<-compute_shrinkage_metrics(basis.DT)
  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  saveRDS(list(basis=pc.emp,shrink=shrink.DT),file=args$cache_file)
}

## perhaps add options for different manifests etc.
option_list = list(
        make_option(c("-s", "--support_dir"), type="character",default=DEFAULT.SUPPORT.DIR,
                    help="Location of support files", metavar="character"),
        make_option(c("-t", "--trait"), type="character",
                    help="trait to simulate, must be in manifest", metavar="character"),
        make_option(c("-n", "--n_sims"), type="numeric",default = DEFAULT.N.SIMS,
                    help="number of simulations to compute", metavar="character"),
        make_option(c("-j", "--snp_stats_dir"), type="character", default=DEFAULT.SNPSTATS.DIR,
                    help="snpStats directory", metavar="character"),
        make_option(c("-p", "--prefix"), type="character", default='',
                    help="output directory", metavar="character"),
        make_option(c("-g", "--gwas_data_dir"), type="character", default=DEFAULT.GWAS.DATA.DIR,
                    help="location of OR aligned GWAS source files", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                    help="output directory", metavar="character"),
        make_option(c("-c", "--cache_file"), type="character", default=NULL,
                    help="file to cache basis in", metavar="character")
        )


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  print(args)
  if(args$prefix=='')
    args$prefix <- basename(tempfile(pattern=sprintf("%s_",args$trait)))
}else{
  message("IN TESTING MODE ======>!")
  trait <- 'ret_p'
  args <- list(
    support_dir = DEFAULT.SUPPORT.DIR,
    trait = trait,
    snp_stats_dir = DEFAULT.SNPSTATS.DIR,
    prefix = basename(tempfile(pattern=sprintf("%s_",trait))),
    n_sims = 2,
    out_dir = '~/tmp/test_distance/',
    gwas_data_dir = DEFAULT.GWAS.DATA.DIR,
    cache_file = '~/tmp/feb_2018_as_basis_cache.RDS'
  )
  print(args)
}


## build basis  -  cache this as does not change between runs
support.dir <- args$support_dir
if(!file.exists(args$cache_file)){
  cacheBasis(args)
  message(sprintf('cached basis in %s',args$cache_file))
  stop()
}else{
  cache.obj <- readRDS(args$cache_file)
}

## simulate
support.dir <- args$support_dir
ref_af_file<-file.path(support.dir,DEFAULT.REF.AF.FILE)
ld_file<-file.path(support.dir,DEFAULT.LD.BLOCK.FILE)
m_file<-file.path(support.dir,DEFAULT.MANIFEST.FILE)
study1.DT<-get_gwas_data(m_file,ref_af_file,args$gwas_data_dir,FALSE,args$trait)
study1.DT <- study1.DT[cache.obj$shrink[,.(pid,emp_maf_se)]]
## simulate under the null
## need tp compute empirical SE for study
study1.DT[,emp_se := 1/sqrt(2) * 1/sqrt(n) * emp_maf_se]
study1.DT[,or:=1]
## need chr and perhaps position
study1.DT[,c('chr','position'):=tstrsplit(pid,':')]
message(sprintf("Simulate %d studies",args$n_sims))
s1.sim <- simulate_study(study1.DT,args$snp_stats_dir,shrink_beta=FALSE,n_sims=args$n_sims,quiet=TRUE)

## compute the projection on 10 at a time to prevent huge matrix generation eating all memory
idx <- split(1:nrow(s1.sim),s1.sim$trait)
pred.emp<-do.call('rbind',lapply(split(idx, ceiling(seq_along(idx)/10)),function(i){
  tmp.idx <- unlist(i)
  tmp.sim <- s1.sim[tmp.idx,]
  setkey(tmp.sim,pid)
  sim.mat.emp<-create_ds_matrix(tmp.sim,cache.obj$shrink,'emp')
  predict(cache.obj$basis,newdata=sim.mat.emp)
}))
ofile <- paste(file.path(args$out_dir,args$prefix),'RDS',sep='.')
saveRDS(pred.emp,file=ofile)
message(sprintf("Saved projections in %s",ofile))
