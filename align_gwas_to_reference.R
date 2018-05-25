## a script to align a GWAS to a reference
## everntually move functions into cupcake
library(data.table)


TEST <- TRUE

option_list = list(
  make_option(c("-g", "--gwas_file"), type="character",default='',
              help="Path to GWAS file", metavar="character"),
  make_option(c("-r", "--reference_file"), type="character",default='',
              help="Path to reference file", metavar="character")
  make_option(c("-o", "--output_file"), type="character",default='',
              help="Path to output file", metavar="character")
)

if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  print(args)
  if(!file.exists(args$gwas_file))
    stop(sprintf("Cannot locate gwas file %s",args$gwas_file))
  if(file.exists(args$reference_file))
      stop(sprintf("Cannot locate reference file %s",args$=reference_file))
  if(file.exists(args$output_file))
    stop(sprintf("Output file already exists %s",args$output_file))
}else{
  message("IN TESTING MODE ======>!")
  trait <- 'jia_cc'
  args <- list(
    target_dir = '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned',
    output_dir = '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned/split/',
    chunk_size = 50
  )
}
