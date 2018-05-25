## It makes sense to compute individual level posterior OR levels for all individuals across a chromosome, however when it comes to projection the matrix to be computed is huge and has a significant but unnecessary burden. This script takes a dir of projected lor objects and an output dir (by default a 'split' dir in the target dir).

library(optparse)
library(data.table)

TEST <- TRUE
DEFAULT_TARGET_DIRNAME <- 'split'

option_list = list(
  make_option(c("-t", "--target_dir"), type="character",default='',
              help="Location of support files", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character",
              help="Location for chunked files to be written", metavar="character"),
  make_option(c("-c", "--chunk_size"), type="numeric",default=50,
              help="Number of samples to use for each chunk", metavar="character")
)


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  print(args)
  if(!file.exists(args$target_dir))
    stop(sprintf("Cannot locate target_dir %s",args$target_dir))
  if(args$output_dir=='')
    args$output_dir <- file.path(args$target_dir,DEFAULT_TARGET_DIRNAME)
}else{
  message("IN TESTING MODE ======>!")
  trait <- 'jia_cc'
  args <- list(
    target_dir = '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned_fixed',
    output_dir = '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned_fixed/split/',
    chunk_size = 50
  )
}

if(!file.exists(args$output_dir)){
  message(sprintf("Creating %s",args$out_dir))
  dir.create(args$output_dir)
}else{
  message(sprintf("Output dir %s",args$output_dir))
}


## next loop through all files in this dir creating split files in the target directory

files <- list.files(path=args$target_dir,pattern="*.RDS",full.names=TRUE)
for(f in files){
  message(sprintf("Processing %s",f))
  chrname <- gsub("([^.]+)\\.RDS","\\1",basename(f))
  tobj <- readRDS(f)
  ## annoyingly when created these objects use factors which makes things slow...
  ## need to fix this when we next recreate
  if(!exists("chunks")){
    all.samples <- as.character(tobj$samples$ID_1)
    chunks <- split(all.samples,ceiling(seq_along(all.samples)/args$chunk_size))
  }
  for(i in seq_along(chunks)){
    sid <- chunks[[i]]
    idx <- which(tobj$samples$ID_1 %in% sid)
    samples <- tobj$samples[idx,]
    samples$ID_1 <- as.character(samples$ID_1)
    samples$ID_2 <- as.character(samples$ID_2)
    samples <- data.table(samples)
    nobj<-list(
      snps = copy(tobj$snps),
      proj_lor = tobj$proj.lor[,idx],
      samples = samples
    )
    ofile <- file.path(args$output_dir,sprintf("%d_%s_tmp.RDS",i,chrname))
    saveRDS(nobj,file=ofile)
    message(sprintf("Processed %s",ofile))
  }
}

## next we reassemble so we have complete genome chunks across all individuals that we can analyse

for(j in seq_along(chunks)){
  pat <- sprintf("^%d_chr.*_tmp.RDS",j)
  cfiles <- list.files(path=args$output_dir,pattern=pat,full.names=TRUE)
  objs <- lapply(cfiles,readRDS)
  nobj<-list(
    snps=rbindlist(lapply(objs,'[[','snps')),
    proj_lor=do.call('rbind',lapply(objs,'[[','proj_lor')),
    samples=objs[[1]]$samples
  )
  ofile <- file.path(args$output_dir,sprintf("%d_all_chr.RDS",j))
  saveRDS(nobj,file=ofile)
  message(sprintf("Processed %s",ofile))
  for(tfiles in cfiles)
    file.remove(tfiles)
}



### just test if this code works


library(cupcake)
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


input.dir <- '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned/split'
in.files <- list.files(path=input.dir,pattern='*.RDS',full.names=TRUE)

library(parallel)
all.preds <- mclapply(in.files,function(f){
  message(f)
  obj <- readRDS(f)
  snps <- obj$snps
  or <- data.table(exp(obj$proj_lor))
  or[,snp:=snps$pid]
  samples <- obj$samples
  setnames(or,c(as.character(samples$ID_1),'pid'))
  DT <- melt(or,id.vars='pid')
  setnames(DT,c('pid','sample','or'))
  setcolorder(DT,c('sample','pid','or'))
  setkey(DT,pid)
  setnames(DT,'sample','trait')
  mat.emp <- create_ds_matrix(DT,shrink.DT,'emp')
  if(!identical(colnames(mat.emp),colnames(basis.mat.emp)))
    stop("Something wrong basis and projection matrix don't match")
  list(sample=obj$samples,pred=predict(pc.emp,newdata=mat.emp))
},mc.cores=8)


all.proj<-data.table(melt(do.call('rbind',lapply(all.preds,'[[','pred'))))
setnames(all.proj,c('sample','pc','value'))
all.samples<-rbindlist(lapply(all.preds,'[[','sample'))
setkey(all.proj,sample)
setkey(all.samples,ID_1)
all.proj <- all.proj[all.samples[,.(ID_1,alt_ilar_code)]]
setnames(all.proj,'alt_ilar_code','trait')
setcolorder(all.proj,c('trait','pc','value','sample'))


jia_traits<-fread(m_file)[grep('jia_',trait),]$trait
#bb_traits<-bb_traits[bb_traits != 'jia_cc']
jia.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,jia_traits)
jia.mat.emp <- create_ds_matrix(jia.DT,shrink.DT,'emp')
jia.pred <- data.table(melt(predict(pc.emp,newdata=jia.mat.emp)))
setnames(jia.pred,c('trait','pc','value'))
jia.pred[,trait:=gsub("jia\\_","",trait)]
comb<-rbind(jia.pred,all.proj,fill=TRUE)
comb[,class:=ifelse(is.na(sample),'summary','individual')]

## really want to compare control
control <- data.table(melt(pc.emp$x))[Var1=='control',]
setnames(control,c('trait','pc','value'))
control[,class:='control']
control<-rbindlist(lapply(unique(comb$trait),function(tra){
  tmp<-copy(control)
  tmp[,trait:=tra]
}))

final <- rbind(control,comb,fill=TRUE)
final[is.na(sample),sample:=paste(trait,class,sep='_')]

library(ggplot2)
 ggplot(final[!trait %in% c('cc','missing','UnA'),],aes(x=pc,y=as.numeric(value),alpha=class!='individual',colour=class,group=sample)) + geom_point() + geom_line() + facet_grid(trait~.) + theme_bw() + scale_alpha_discrete(guide=FALSE)


pc <- dcast(final,trait+class+sample~pc)[!trait %in% c('cc','missing','UnA'),]

 ggplot(pc,aes(x=PC1,y=PC2,alpha=class!='individual',colour=class)) + geom_point()  + facet_grid(trait~.) + the

 ggplot(pc[trait=='ERA',],aes(x=PC1,y=PC2,alpha=class!='individual',colour=class)) + geom_point() + theme_bw()

 ## perhaps split is screwing us up ?

 out_dir <- '/scratch/ob219/as_basis/jia_ind_analysis/ind_proj_aligned'
 dat<-readRDS(file.path(out_dir,"chr22.RDS"))
 sman<-dat$snps
 sman[,mean.proj.lor:=rowMeans(dat$proj.lor)]

 ## next do the same with our splits but for same chr22

 all.sums <- mclapply(in.files,function(f){
   message(f)
   obj <- readRDS(f)
   snps <- obj$snps
   chr.idx<-which(snps$chr=='22')
   lor <- data.table(obj$proj_lor)
   lor <- lor[chr.idx,]
},mc.cores=8)
all.sums<-do.call(cbind,all.sums)
obj <- readRDS(in.files[1])
chr.idx<-which(obj$snps$chr=='22')
snps <- obj$snps[chr.idx,]
snps[,mean.proj.lor:=rowMeans(all.sums)]
snps[flip==TRUE,beta:=-beta]
ggplot(snps,aes(x=beta,y=mean.proj.lor)) + geom_point() + geom_smooth(method="lm") +
geom_abline(intercept=0,slope=1,color='red',lty=2) + theme_bw() + xlab("Case/Control Beta") +
ylab("Mean(Beta.proj)")

#seems OK what about the OR from the 'aligned GWAS'

cc <- jia.DT[trait=='jia_cc' & chr=='22',.(pid,basis.beta=log(or))]
setkey(cc,pid)
setkey(snps,pid)
snps<-snps[cc]
ggplot(snps,aes(x=basis.beta,y=mean.proj.lor,color=flip)) + geom_point() + geom_smooth(method="lm") +
geom_abline(intercept=0,slope=1,color='red',lty=2) + theme_bw() + xlab("Case/Control Beta") +
ylab("Mean(Beta.proj)")

## plot summary OR against each other

ggplot(snps,aes(x=basis.beta,y=beta,color=flip)) + geom_point()


## OK have messed up - need to fix the underlying files so that OR are correct wrt to the basis

obj <- readRDS(in.files[1])
snps <- obj$snps
snps[flip==TRUE,beta:=-beta]
cc <- jia.DT[trait=='jia_cc',.(pid,basis.beta=log(or))]
setkey(cc,pid)
setkey(snps,pid)
snps<-snps[cc]
ggplot(snps,aes(x=basis.beta,y=beta,color=flip)) + geom_point()
## ok should be x==y
snps[,poss.flip:=sign(beta)!=sign(basis.beta)]
ggplot(snps,aes(x=basis.beta,y=beta,color=poss.flip)) + geom_point()

## double check that we did not inadvertantly mess up the alleles

## so we want to compute the mean OR over all individuals for a particular subtype

all.sums <- mclapply(in.files,function(f){
  message(f)
  obj <- readRDS(f)
  snps <- obj$snps
  or <- data.table(exp(obj$proj_lor))
  or[,snp:=snps$pid]
  samples <- obj$samples
  setnames(or,c(as.character(samples$ID_1),'pid'))
  DT <- melt(or,id.vars='pid')
  setkey(DT,variable)
  setkey(samples,ID_1)
  DT <- DT[samples[,.(ID_1,alt_ilar_code)]]
  DT[,list(sor=sum(value),n.samples=.N),by=c('pid','alt_ilar_code')]
},mc.cores=8)

all.snps<-readRDS(in.files[1])$snps

all.sums <- rbindlist(all.sums)
final.sums<-all.sums[,list(mean_or=sum(sor)/sum(n.samples)),by=c('pid','alt_ilar_code')]
## next merge in summary or
setkeyv(final.sums,c('pid','alt_ilar_code'))
jia.DT[,alt_ilar_code:=gsub("jia\\_","",trait)]
setkeyv(final.sums,c('pid','alt_ilar_code'))
setkeyv(jia.DT,c('pid','alt_ilar_code'))
lor.plot <- final.sums[jia.DT[,.(pid,alt_ilar_code,or=or)]][!is.na(mean_or),]
lor.plot[,flip:=pid %in% all.snps[flip==TRUE,]$pid]
foo<-lor.plot[alt_ilar_code=='ERA' & or<3,]
foo<-foo[sample.int(nrow(foo),1000,replace=FALSE),]
foo<-foo[,or:=1/or]
ggplot(foo,aes(x=or,y=mean_or,color=flip)) + geom_point()

ggplot(lor.plot,aes(x=lor,y=mean_lor)) + geom_hex() + facet_grid(alt_ilar_code~.)


## scree plots for each


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

dat <- get(load('/scratch/wallace/JIA-2017-data/annotsnpstats-22.RData'))
sano <- data.table(samples(dat))[phenotype==1,.(ID_1,recodedilarcode,alt_ilar_code)]
phe.lu <- split(sano$alt_ilar_code,as.character(sano$ID_1))
emp[11:nrow(emp),compare:=unlist(phe.lu[trait])]
emp[11:nrow(emp),trait:='']
emp[1:9,compare:='none']


ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')
emp[!emp$compare %in% c('control','ERA'),compare:='other.subtype']


jia.only<-subset(emp,!compare %in% c('control','none','missing'))
t.test(jia.only[jia.only$compare=='ERA',]$PC2,jia.only[jia.only$compare!='ERA',]$PC2)
