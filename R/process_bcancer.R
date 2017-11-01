# ## process breast cancer data.
#
# bc<-fread('oncoarray_bcac_public_release_oct17.txt')
#
# dtypes <- rev(c('P1df','se','beta','eaf_controls','P1df_Wald','P1df_LRT','r2'))
#
# ## we want to split into different GWAS
#
# stub<-names(bc)[1:6]
# dnames<-names(bc)[7:length(names(bc))]
# ## to collect use dtypes
# for(d in dtypes){
# #lapply(dtypes,function(d){
#   dnames <- gsub(paste0('_',d),'',dnames)
# }
#
# out.dir <- '/home/ob219/scratch/as_basis/gwas_stats/raw/bcac/'
#
# by.an <- split(1:length(dnames),dnames)
#
# lapply(seq_along(by.an),function(i){
#   ci<-by.an[[i]]
#   fname<-file.path(out.dir,sprintf("%s.RDS",names(by.an)[i]))
#   saveRDS(bc[,c(1:6,ci),with=FALSE],file=fname)
# })
#
# ## next we do the above but just for basis SNPs
#
# bc[,pid:=paste(chr,position_b37,sep=':')]

## get a list of basis SNPs


#lapply(seq_along(by.an),function(i){
#  ci<-by.an[[i]]
#  fname<-file.path(out.dir,sprintf("filtered_%s.RDS",names(by.an)[i]))
#  saveRDS(bf[,c(1:6,ci),with=FALSE],file=fname)
#})

## align alleles with the basis for each trait and format for cupcake

## get the format
library(data.table)
dat.file <- '/home/ob219/scratch/as_basis/gwas_stats/raw/bcac/oncoarray_bcac_public_release_oct17.txt'
bc <- fread(dat.file,nrows=1L)
bc<-fread(dat.file,select=names(bc)[1:10])
bc[,pid:=paste(chr,position_b37,sep=':')]
## what is N ?

man<-data.table(study=c('BCAC','ICOGS','GWAS'),cases=c(61282,46785,14910),controls=c(45494,42892,17588))

n.cases<-sum(man$cases)
n.controls<-sum(man$controls)

support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snps_with_alleles.tab')
bsnps<-fread(ref_af_file)
bsnps[,pid:=paste(chr,position,sep=':')]
setkey(bsnps,pid)
bf<-bc[pid %in% bsnps$pid,]
bf<-bf[!duplicated(pid),]
setkey(bf,pid)

bfo<-bf[,.(id=phase3_1kg_id,chr,position=position_b37,a1=a0,a2=a1,or=exp(as.numeric(bcac_onco_icogs_gwas_beta)),p.val=bcac_onco_icogs_gwas_P1df)]




library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)


res<-align_alleles(bfo,bsnps,check=TRUE)
## or are flipped thus we need 1/or !!
res<-res[,.(id,chr,position,p.val=sprintf("%.4f",p.val),or=sprintf("%.4f",1/or))]
#res<-res[,.(id,chr,position,p.val=sprintf("%.4f",p.val),or=sprintf("%.4f",or))]
options(scipen=999)
write.table(res,file='/home/ob219/scratch/as_basis/gwas_stats/input_files/breast_cancer.tab',sep="\t",row.names=FALSE,quote=FALSE)
options(scipen=0)
