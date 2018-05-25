library(data.table)
library(magrittr)
library(cupcake)
library(ggplot2)
library(annotSnpStats)
# aggregate chromosomes across split sample files

in.dir <- '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned'
files <- list.files(path=in.dir,pattern='*.RDS',full.names=TRUE)
by.run <- split(files,strsplit(basename(files),'_') %>% sapply(.,'[[',2))

## compute the basis

support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snps.tab')
ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-file.path(support.dir,'as_basis_manifest_with_jia_cc.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

library(parallel)
all.pred <- mclapply(seq_along(by.run),function(j){
  message(sprintf("Processing %s",names(by.run)[j]))
  fs <- by.run[[j]]
  DT <- rbindlist(lapply(fs,readRDS))
  summ <- DT[,list(total=sum(or),n=.N),by=pid]
  ## for speed do 10 samples at a time
  setnames(DT,'sample','trait')
  setkey(DT,trait)
  samples<-unique(DT$trait)
  sc <- split(samples, ceiling(seq_along(samples)/10))
  do.call('rbind',(lapply(seq_along(sc),function(i){
    #message(i)
    s <- sc[[i]]
    tmp.DT <- subset(DT,trait %in% s)
    ## test what happens if we swapped alleles incorrectly
    tmp.DT[flipped=-1,or:=1/or]
    ## assume that SNPs are in the same order for time being
    mat.emp <- create_ds_matrix(tmp.DT,shrink.DT,'emp')
    if(!identical(colnames(mat.emp),colnames(basis.mat.emp)))
      stop("Something wrong basis and projection matrix don't match")
    pred.emp <- predict(pc.emp,newdata=mat.emp)
    return(pred.emp)
    #emp<-rbind(pc.emp$x,pred.emp)
  })))
},mc.cores=8)

all.proj<-do.call('rbind',all.pred)

bb_traits<-fread(m_file)[grep('jia_',trait),]$trait
#bb_traits<-bb_traits[bb_traits != 'jia_cc']
bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,bb_traits)

emp<-rbind(pc.emp$x,all.proj)
ml<-list(
  CD = 'bb_CD',
  CEL = 'bb_CEL',
  MS = 'bb_MS',
  RA = 'bb_RA',
  SLE = 'bb_SLE',
  T1D = 'bb_T1D',
  UC = 'bb_UC'
)
g <- function(M){
    M <- cbind(as.data.table(M),trait=rownames(M))
    M$compare<-"none"
    for(i in seq_along(ml)) {
        M[trait %in% c(names(ml)[i], ml[i]), compare:=names(ml)[i]]
    }
    M[trait=="control",compare:="control"]
    M
}
emp<-g(emp)





## need to get the subtype data

dat <- get(load('/scratch/wallace/JIA-2017-data/annotsnpstats-22.RData'))
sano <- data.table(samples(dat))[phenotype==1,.(ID_1,recodedilarcode,alt_ilar_code)]
phe.lu <- split(sano$alt_ilar_code,as.character(sano$ID_1))
emp[11:nrow(emp),compare:=unlist(phe.lu[trait])]
emp[11:nrow(emp),trait:='']
emp[1:9,compare:='none']


ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')
emp[!emp$compare %in% c('control','ERA'),compare:='other.subtype']

ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')
jia.only<-subset(emp,!compare %in% c('control','none','missing'))

t.test(jia.only[jia.only$compare=='ERA',]$PC2,jia.only[jia.only$compare!='ERA',]$PC2)

snp.summ <- mclapply(seq_along(by.run),function(j){
  message(sprintf("Processing %s",names(by.run)[j]))
  fs <- by.run[[j]]
  tf <- rbindlist(lapply(fs,function(f){
    message(sprintf("Processing %s",f))
    DT <- readRDS(f)
    DT[,ilar:=unlist(phe.lu[sample])]
    DT[,list(total=sum(or),n=.N),by=c('pid','ilar')]
  }))
},mc.cores=8)

ilar.mean<-rbindlist(snp.summ)[,list(mean=sum(total)/sum(n)),by=c('pid','ilar')]
setkeyv(ilar.mean,c('pid','ilar'))

bb_traits<-fread(m_file)[grep('jia_',trait),]$trait
#bb_traits<-bb_traits[bb_traits != 'jia_cc']
bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,bb_traits)
bb.DT[,ilar:=sub('jia_','',trait)]
setkeyv(bb.DT,c('pid','ilar'))

pm<-ilar.mean[bb.DT[,.(pid,ilar,or)]][!is.na(mean),]

library(ggplot2)

## remove outliers

rm.pid <- unique(pm[pm$or>2 | pm$or<0.5,]$pid)
pm <- pm[!pid %in% rm.pid,]

ggplot(pm,aes(x=log(or),y=log(mean))) + geom_hex() + facet_grid(.~ilar) + geom_smooth()

## checkout why alleles appear switched


switchy <- subset(pm,ilar=='PsA' & (round(log(mean),digits=3)==-round(log(or),digits=3)))

ggplot(switchy,aes(x=log(or),y=log(mean))) + geom_point()

## looks like we have them

switchy$pid

bsnps<-fread(ref_af_file)[,pid:=paste(chr,position,sep=':')]
## no MAF or location pattern
tp <- bsnps[pid %in% switchy$pid,][order(chr,position),csum:=cumsum(as.numeric(position))]
ggplot(tp,aes(x=csum,y=maf))+ geom_point()

## get the annotation for all SNPs.

if(FALSE){
  jfil <- list.files(path="/scratch/wallace/JIA-2017-data",pattern="annotsnpstats*",full.names=TRUE)
jia.snps <- mclapply(jfil,function(file){
  message("Processing %s",file)
  D<-get(load(file))
  data.table(snps(D))
},mc.cores=8)
all.snps<-rbindlist(jia.snps)
saveRDS(all.snps,'/scratch/wallace/JIA-2017-data/snpAnnotations.RDS')
}
else{
  all.snps <- readRDS('/scratch/wallace/JIA-2017-data/snpAnnotations.RDS')
}

foo <- subset(all.snps,pid %in% switchy$pid)

foo[,flip:=allele.1==ALT & allele.2==REF] ## could be this ?

## get a list of all those that might need to be flipped

jia.BT.DT <- all.snps[pid %in% bsnps$pid,][,flip:=allele.1==ALT & allele.2==REF]
setkey(jia.BT.DT,pid)
setkey(pm,pid)
pm <- pm[jia.BT.DT[,.(pid,flip)]]
pm[,mean.cor:=mean]
pm[flip==TRUE,mean.cor:=1/mean]

ggplot(pm,aes(x=log(or),y=log(mean.cor))) + geom_hex() + facet_grid(.~ilar)



pm[,mean.cor:=1/or]
