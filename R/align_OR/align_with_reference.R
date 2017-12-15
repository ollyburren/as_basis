## code to align a dataset with a reference


library(optparse)
library(data.table)
library(magrittr)

TEST<-FALSE
MISSING_FILE_DEFAULT<-'/scratch/ob219/as_basis/support_tab/additional_missing_from_basis.txt'
option_list = list(
        make_option(c("-i", "--input"), type="character", default=NULL,
              help="input tab file to process", metavar="character"),
        make_option(c("-o", "--output"), type="character", default=NULL,
              help="output tab file to process", metavar="character"),
        make_option(c("-r", "--ref"), type="character", default=NULL,
              help="reference file to use", metavar="character"),
        make_option(c("-a", "--effect_allele"), type="character", default=NULL,
              help="Allele which OR is with respect to",metavar="character"),
        make_option(c("-m", "--missing_file"), type="character", default=MISSING_FILE_DEFAULT,
              help="Were to write SNPs missing from existing reference file",metavar="character")
        )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}else{
  args <- list(
      ref='/scratch/ob219/as_basis/support_tab/as_basis_snp_support.tab',
      output='~/tmp/test_align.tab',
      input='/scratch/ob219/as_basis/gwas_stats/processed//bb:20002_1065:self_reported_hypertension.tab',
      effect_allele='a2',
      missing_file=MISSING_FILE_DEFAULT
  )
}

message(sprintf("##############\nProcessing %s\n##############",args$input))
ref <- fread(args$ref)
setkey(ref,pid)
setnames(ref,c('pid','ref_a1','ref_a2','ref_a1.af','ld.block'))
input <- fread(args$input)[nchar(a1)==1 & nchar(a2)==1,]
setkey(input,pid)
out <- input[ref][!is.na(ref_a1.af),]
## convert so OR is always wrt a2
if(args$effect_allele=='a1')
  out[,or:=1/or]

comp<-function(cv){
  tmp<-cv
  a<-c('A','G','C','T')
  b<-c('T','C','G','A')
  idx<-split(1:length(cv),cv)
  for(i in seq_along(a)){
    cv[idx[[a[i]]]]<-b[i]
  }
  cv
}

## now attempt to align

out[,flag:='unprocessed']
out[a1==ref_a1 & a2==ref_a2,flag:='match']
out[a1==ref_a2 & a2==ref_a1,flag:='flip']
out[a1==comp(ref_a1) & a2==comp(ref_a2),flag:='match_rc']
out[a1==comp(ref_a2) & a2==comp(ref_a1),flag:='flip_rc']

if(any(out$flag=='unprocessed')){
  probs <- out[flag=='unprocessed',]$pid
  if(length(unique(out[flag=='unprocessed',]$pid))==length(unique(out[flag!='unprocessed',]$pid))){
    print(out[flag=='unprocessed',])
    stop("Issue with matching alleles to reference")
  }
}

## for flips we need to flip or
if(any(out$flag=='unprocessed')){
  # we should store these in the missing file
  missing <- unique(out[flag=='unprocessed',]$pid)
  if(file.exists(args$missing_file)){
    missing <- unique(union(scan(args$missing_file,character()),missing))
  }
  write(missing,file=args$missing_file)
}
out <- out[flag!='unprocessed',]
if(any(out$flag %in% c('flip','flip_rc')))
  out[flag %in% c('flip','flip_rc'),c('ori.or','or'):=list(or,1/or)]
message(sprintf("writing %d rows to %s",nrow(out),args$output))
out <- out[,.(pid,ref_a1,ref_a2,or=signif(or,digits=5),p.value=signif(p.value,digits=5))]
setnames(out,c('pid','a1','a2','or','p.value'))
write.table(out,file=args$output,quote=FALSE,row.names=FALSE,sep='\t')

## this to process the whole manifest

if(FALSE){
  library(data.table)
  rscript <- '/home/ob219/git/as_basis/R/align_OR/align_with_reference.R'
  process.dir<-'/scratch/ob219/as_basis/gwas_stats/processed/'
  manifest <- fread('/scratch/ob219/as_basis/support_tab/as_manifest_december.tab')
  out.dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files/'
  ref<-'/scratch/ob219/as_basis/support_tab/as_basis_snp_support.tab'
  all.cmds <- lapply(1:nrow(manifest),function(i){
    template <- "Rscript --vanilla %s -r %s -i %s -a %s -o %s"
    ifile <- file.path(process.dir,sprintf("%s.tab",manifest$disease[i]))
    ofile<-file.path(out.dir,basename(ifile))
    ef <- manifest$effect_allele[i]
    sprintf(template,rscript,ref,ifile,ef,ofile)
  })
  write(do.call('c',all.cmds),file='/home/ob219/git/as_basis/sh/align_alleles.sh')

}
