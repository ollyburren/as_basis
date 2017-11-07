library(data.table)
library(magrittr)
library(cupcake)
# aggregate chromosomes across split sample files

in.dir <- '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_basis'
files <- list.files(path=in.dir,pattern='*.RDS',full.names=TRUE)
by.run <- split(files,strsplit(basename(files),'_') %>% sapply(.,'[[',2))

## compute the basis

support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snps.tab')
ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-file.path(support.dir,'as_basis_manifest.tab')
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
all.pred <- mclapply(seq_along(head(by.run,n=1)),function(j){
  message(sprintf("Processing %s",names(by.run)[j]))
  fs <- by.run[[j]]
  DT <- rbindlist(lapply(fs,readRDS))
  ## for speed do 10 samples at a time
  setnames(DT,'sample','trait')
  setkey(DT,trait)
  samples<-unique(DT$trait)
  sc <- split(samples, ceiling(seq_along(samples)/10))
  do.call('rbind',(lapply(seq_along(sc),function(i){
    #message(i)
    s <- sc[[i]]
    tmp.DT <- subset(DT,trait %in% s)
    setkey(tmp.DT,pid)

    ## test what happens if we swapped alleles incorrectly
    tmp.DT[,or:=1/or]

    ## assume that SNPs are in the same order for time being
    mat.emp <- create_ds_matrix(tmp.DT,shrink.DT,'emp')
    if(!identical(colnames(mat.emp),colnames(basis.mat.emp)))
      stop("Something wrong basis and projection matrix don't match")
    pred.emp <- predict(pc.emp,newdata=mat.emp)
    return(pred.emp)
    #emp<-rbind(pc.emp$x,pred.emp)
  })))
},mc.cores=4)

all.proj<-do.call('rbind',all.pred)

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

bDT <- tmp.DT
sDT <- shrink.DT
method <- 'emp'

create_ds_matrix <- function(bDT,sDT,method=c('emp','est')){
  if(missing(method)){
    method='emp'
  }
  message(sprintf("Using %s",method))
  vmethod = sprintf("%s_shrinkage",method)
  stmp<-sDT[,c('pid',vmethod),with=FALSE]
  tmp<-bDT[stmp]
  tmp$metric <- tmp[[vmethod]] * log(tmp$or)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  return(as.matrix(B[,-1]) %>% t())
}


## need to get the subtype data
library(annotSnpStats)
dat <- get(load('/scratch/wallace/JIA-2017-data/annotsnpstats-22.RData'))
sano <- data.table(samples(dat))[phenotype==1,.(ID_1,recodedilarcode,alt_ilar_code)]
phe.lu <- split(sano$alt_ilar_code,as.character(sano$ID_1))
emp[11:nrow(emp),compare:=unlist(phe.lu[trait])]
emp[11:nrow(emp),trait:='']
emp[1:9,compare:='none']


ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')

jia.only<-subset(emp,!compare %in% c('control','none','missing'))

t.test(jia.only[jia.only$compare=='ERA',]$PC2,jia.only[jia.only$compare!='ERA',]$PC2)
