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
## code for converting linear regression coefficient ot an OR
#source('/home/ob219/git/as_basis/R/convertBetaToOR.R')

computeOR <- function(n,n1,Sx,Sxy) {
    ## estimated allele freqs in cases and controls
    fe1 <- Sxy/(2*n1)
    fe0 <- (Sx - Sxy)/(2*(n-n1))
    ## estimated odds ratio
    fe1 * (1-fe0) / ( (1-fe1) * fe0 )
}

SElor<-function(n,n1,Sx,Sxy){
    n0<-n-n1
    #fe1 is the af in cases
    c <- Sxy/(2*n1)
    #fe0 is af in the controls
    a <- (Sx - Sxy)/(2*(n-n1))
    b<-1-a
    d<-1-c
    ## normalise
    a<-(a*n0)/n
    b<-(b*n0)/n
    c<-(c*n1)/n
    d<-(d*n1)/n
    ## estimated odds ratio bc/ad
    sqrt(1/2) * sqrt(1/n) * sqrt(1/a + 1/b + 1/c + 1/d)
}


data.dir<-'/home/ob219/scratch/as_basis/bb/summary_stats/'
out.dir<-'/home/ob219/scratch/as_basis/gwas_stats/processed/'
bb_phenofile<-'/home/ob219/scratch/as_basis/bb/phenosummary_final_11898_18597.tsv'
all.summ.files<-list.files(path=data.dir,pattern='*.RDS',full.names=TRUE)

## this is a list of the SNPs that are in the basis - we don't care about SNPs that are not in this list

#basis.snps<-readRDS('/home/ob219/scratch/as_basis/basis_snp_positions.RDS')
bb_var_info<-readRDS('/home/ob219/scratch/as_basis/bb/variants_prefiltered.RDS')
## read in phenotype file
pheno<-fread(bb_phenofile)
## should probably do this on the queue as these are very large files

f<-args$file
## get field code
field.code <- unlist(strsplit(basename(f),':'))[1]
message(sprintf("Loading %s",f))
stats<-readRDS(f)
message(sprintf("Finshed Loading %s",f))
#stats[,lu:=sub('(.*:.*):.*','\\1',variant)]
n.cases<-subset(pheno,Field.code==field.code)$N.cases
#lu<-substr(stats$variant,1,nchar(stats$variant)-4)
keep<-which(stats$variant %in% bb_var_info$variant)

## some SNPs are missing - just need to make sure that these are removed when computing the basis.
#miss<-basis.snps[-which(basis.snps %in% lu)]

stats.filter<-stats[keep,]
stats.filter[,c('chr','position','A1','A2'):=tstrsplit(variant, ":", fixed=TRUE)]

## add MAF so we can compute an OR
bb_var_info<-bb_var_info[,.(variant,AF)]
setkey(bb_var_info,variant)
setkey(stats.filter,variant)
stats.filter<-bb_var_info[stats.filter]

# this readme says that the effect allele is the same as ref
#https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/README.md#summary-stat-

## these are based on a linear regression so we need to alter so we can get an OR
stats.filter[,maf:=ifelse(AF>0.5,1-AF,AF)]
stats.filter[,or:=computeOR(nCompleteSamples,n.cases,AC,ytx)]
stats.filter[,c('theta','se.theta'):=list(log(or),SElor(nCompleteSamples,n.cases,AC,ytx))]
stats.filter[,c('theta.pval','theta.Z','n0','n1'):=list(2*(pnorm(abs(theta/se.theta),lower.tail = FALSE)),theta/se.theta,nCompleteSamples-n.cases,n.cases)]

## create the classic output that we need to include into the basis

#id chr position risk.allele other.allele or p.val

stats.filter[,c('risk.allele','other.allele','id'):=list(A1,A2,paste(chr,position,sep=':'))]
## find those that we need to flip
stats.filter[which(stats.filter$or<1),c('or','risk.allele','other.allele'):=list(1/or,A2,A1)]
## output

trait<-paste(gsub('[.]+','_',sub('.*:Non\\.cancer\\.illness\\.code\\.\\.(.*)\\.RDS','\\1',basename(f))),'out',sep='.')
setnames(stats.filter,c('theta.pval','rsid'),c('p.val','id'))
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

if(FALSE){
  QQ=function(x,l=0.99,add=FALSE,minx=TRUE,...) {
  if (max(abs(x),na.rm=T)<1.1) x=-log10(x)
  n=length(x); q=-log10((n:1)/(n+1)); x=sort(x)
  n1=round(l*n)
  if (minx) {
    if (add) points(c(0,q[n1:n]),c(0,x[n1:n]),...) else plot(c(0,q[n1:n]),c(0,x[n1:n]),...)
    lines(c(0,q[n1]),c(0,x[n1]),...)
  } else if (add) points(q,x,...) else plot(q,x,...)
}

QQ(stats.filter$pval,l=0.99,xlab='Expected(-log(P))',ylab='Observed(-log(P))',main='QQ plot linear(black) vs logistic(red)')
QQ(stats.filter$theta.pval,l=0.99,add=TRUE,col='red')
}
