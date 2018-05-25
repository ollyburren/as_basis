library(data.table)
library(magrittr)

library(optparse)

TEST <- FALSE

option_list = list(
        make_option(c("-i", "--in_file"), type="character", default=NULL,
              help="GWAS file to filter", metavar="character"),
        make_option(c("-s", "--support_file"), type="character",
                    help="SNP support file to use for filter", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                    help="Directory to write results", metavar="character")
        )

if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if(!file.exists(args$in_file))
    stop("Can't find source data files")
  if(!file.exists(args$support_file))
    stop("No support file")
}else{
  #for testing
  message("IN TESTING MODE ======>!")
  args <- list(
    in_file='/home/ob219/scratch/as_basis/gwas_stats/processed_new/uc_delaange.tab',
    support_file='/home/ob219/scratch/as_basis/support_tab/as_basis_snp_support_feb_2018.tab',
    out_dir='/home/ob219/scratch/as_basis/gwas_stats/filter_feb_2018/'
  )
}

support <- fread(args$support_file)
gwas <- fread(args$in_file)[pid %in% support$pid,]
#gwas <- fread(IN_FILE)
#gwas.filt <- gwas[pid %in% support$pid,]
fout <- file.path(args$out_dir,basename(args$in_file))
if(nrow(support) == nrow(gwas)){
  write.table(gwas,file=fout,quote=FALSE,row.names=FALSE,sep="\t")
  message(sprintf("Written %s",fout))
}else{
  message(sprintf("Something wrong with %s",args$in_file))
}

#if(FALSE){
#  in.files<-list.files(path='/home/ob219/scratch/as_basis/gwas_stats/processed_new/',pattern="*.tab",full.names=TRUE)
#  sfile <- '/home/ob219/scratch/as_basis/support_tab/as_basis_snp_support_feb_2018.tab'
#  odir <- '/home/ob219/scratch/as_basis/gwas_stats/filter_feb_2018/'
#  all.cmds <- sapply(in.files,function(f){
#      cmd <- sprintf("Rscript --vanilla /home/ob219/git/as_basis/R/align_OR/filter_gwas.R -i %s -s %s -o %s",f,sfile,odir)
#  })
#  write(all.cmds,file="~ob219/git/as_basis/sh/run_gwas_filter.txt")
#}

if(FALSE){
  in.files<-list.files(path='/home/ob219/scratch/as_basis/gwas_stats/filter_feb_2018/unaligned/',pattern="*.tab",full.names=TRUE)
  sfile <- '/home/ob219/scratch/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
  odir <- '/home/ob219/scratch/as_basis/gwas_stats/filter_feb_2018_w_ms/unaligned/'
  all.cmds <- sapply(in.files,function(f){
      cmd <- sprintf("Rscript --vanilla /home/ob219/git/as_basis/R/align_OR/filter_gwas.R -i %s -s %s -o %s",f,sfile,odir)
  })
  write(all.cmds,file="~ob219/git/as_basis/sh/run_gwas_filter.txt")
}
