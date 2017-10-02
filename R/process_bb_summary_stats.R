## batch processing of neale lab summary stats for projection onto the basis.

library(optparse)


option_list = list(
        make_option(c("-f", "--file"), type="character", default=NULL,
              help="input RDS file to convert", metavar="character")
            )

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
if (is.null(args$file)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (output file)", call.=FALSE)
}


library(data.table)

data.dir<-'/home/ob219/scratch/as_basis/bb/summary_stats/'
out.dir<-'/home/ob219/scratch/as_basis/gwas_stats/processed/'

all.summ.files<-list.files(path=data.dir,pattern='*.RDS',full.names=TRUE)

## this is a list of the SNPs that are in the basis - we don't care about SNPs that are not in this list

basis.snps<-readRDS('/home/ob219/scratch/as_basis/basis_snp_positions.RDS')

## should probably do this on the queue as these are very large files

f<-args$file
## get field code
field.code <- unlist(strsplit(basename(f),':'))[1]
message(sprintf("Loading %s",f))
stats<-readRDS(f)
message(sprintf("Finshed Loading %s",f))
stats[,lu:=sub('(.*:.*):.*','\\1',variant)]

lu<-substr(stats$variant,1,nchar(stats$variant)-4)
keep<-which(lu %in% basis.snps)

## some SNPs are missing - just need to make sure that these are removed when computing the basis.
miss<-basis.snps[-which(basis.snps %in% lu)]

stats.filter<-stats[keep,]
# this readme says that the effect allele is the same as ref
#https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/README.md#summary-stat-

stats.filter[,c('chr','position','A1','A2'):=tstrsplit(variant, ":", fixed=TRUE)]

## create the classic output that we need to include into the basis

#id chr position risk.allele other.allele or p.val

stats.filter[,c('or','risk.allele','other.allele','id'):=list(exp(beta),A1,A2,paste(chr,position,sep=':'))]
## find those that we need to flip
stats.filter[which(stats.filter$or<1),c('or','risk.allele','other.allele'):=list(1/or,A2,A1)]
## output

trait<-paste(gsub('[.]+','_',sub('.*:Non\\.cancer\\.illness\\.code\\.\\.(.*)\\.RDS','\\1',basename(f))),'out',sep='.')
setnames(stats.filter,c('pval','rsid'),c('p.val','id'))
assign(eval(trait),stats.filter[,.(id,chr,position,risk.allele,other.allele,or,p.val)])
fout<-paste('bb',field.code,gsub('\\.out','.RData',trait),sep=':')
save(list=trait,file=file.path(out.dir,fout))
message(sprintf("Save %s to %s",trait,file.path(out.dir,fout)))

if(FALSE){
  ## code to generate script to run on HPC
  foo<-do.call('c',lapply(list.files(path=data.dir,pattern='*.RDS',full.names=TRUE),function(f){
    sprintf('Rscript --vanilla /home/ob219/git/as_basis/R/process_bb_summary_stats.R -f %s',f)
  }))
  write(foo,file='/home/ob219/git/as_basis/sh/process_bb.txt')
}
