## Code for computing summary statistics for JIA subtypes.

## parse options first in case there is a problem otherwise slow.
library(optparse)
DEBUG=FALSE
option_list = list(
        make_option(c("-o", "--out"), type="character", default=NULL,
              help="output dir", metavar="character"),
        make_option(c("-p", "--phenotype"), type="character", default='cc',
                    help="Phenotype to use, default cc is just JIA vs Controls", metavar="character"),
        make_option(c("-i", "--input"), type="character", default=NULL,
                                help="Input file (must be of type annotSnpStats)", metavar="character")
        )


## valid phenotypes

PHE<-c('EO','ERA','PO','PsA','RFneg','RFpos','sys','cc')

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
if (is.null(args$out)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (output file)", call.=FALSE)
}

if(!args$phenotype %in% PHE){
  print_help(opt_parser)
	stop(sprintf("%s does not appear to be a vaild phenotype",args$phenotype), call.=FALSE)
}

## phenotypes that are available are

#args<-list(
#  input='/scratch/wallace/JIA-2017-data/annotsnpstats-22.RData',
#  out='/scratch/wallace/JIA-2017-data/summary_stats_OB/',
#  phenotype='sys'
#  )

## for chromosome 1 the largest there are 562166 snps it takes approx 5 seconds to fit glm for 100

library(annotSnpStats)
library(data.table)
## get the data
message(sprintf("Getting data from %s",args$input))
as<-get(load(args$input))
message(sprintf("Finished loading data from %s",args$input))
## create our own column for phenotype
message(sprintf("Filter phenotype %s",args$phenotype))
if(args$phenotype != 'cc'){
  ## phenotype is a bit different as need to remove all other cases as these are not controls
  idx<-which(samples(as)$phenotype==0 | samples(as)$alt_ilar_code==args$phenotype)
  as<-as[idx,]
}


#system.time(snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(as),sets=1:100 ,snp.data=as(as,"SnpMatrix")))
#for testing
if(DEBUG){
  message("DEBUG mode on")
  idx<-1:100
}else{
  idx<-1:nrow(snps(as))
}
message(sprintf("Fitting models for %d snps",max(idx)))
test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(as),sets=idx,snp.data=as(as,"SnpMatrix"))
message("Reformatting and computing summary statistics")
## convert this to a dataframe and add annotation information as well as maf
tmp<-data.frame(do.call('rbind',test))
tmp$rsid<-rownames(tmp)
for(n in colnames(tmp)){
  tmp[[n]]<-unlist(tmp[[n]])
}
test<-data.table(tmp)
test[,Y.var:=args$phenotype]

## compute the summary stats
controls.idx<-which(samples(as)$phenotype==0)
cases.idx<-which(samples(as)$phenotype==1)

cas.ss<-col.summary(as[cases.idx,idx])
ctrl.ss<-col.summary(as[controls.idx,idx])

processSummary<-function(df,prefix){
  rsid<-rownames(df)
  DT<-data.table(df)
  DT<-DT[,.(Calls,MAF)]
  setnames(DT,paste(prefix,names(DT),sep='.'))
  DT$rsid<-rsid
  setkey(DT,rsid)
  DT
}

setkey(test,rsid)

test<-test[processSummary(cas.ss,'case')]
test<-test[processSummary(ctrl.ss,'control')]

## add in allele codings

info<-snps(as[,idx])
info<-data.table(info)
## AF here is with respect to ALT allele
info<-info[,.(ID,chromosome,position,allele.1,allele.2,REF,ALT,AF)]
setnames(info,c('rsid','chr','position','a1','a2','ref','alt','af.wrt.alt'))
setkey(info,rsid)
test<-info[test]
## finally compute Z score and p.value

test[,Z:=beta/sqrt(Var.beta)]
test[,p.val:=2*pnorm(abs(Z),lower.tail=FALSE)]
out.file<-sprintf("%d_%s.RDS",unique(test$chr),args$phenotype)
out.file<-file.path(args$out,out.file)
saveRDS(test,file=out.file)
message(sprintf("Wrote file to %s",out.file))
