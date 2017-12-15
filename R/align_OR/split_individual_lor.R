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
ref_af_file<-file.path(support.dir,'as_basis_snp_support.tab')
m_file<-file.path(support.dir,'as_manifest_december.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)


input.dir <- '/home/ob219/scratch/as_basis/jia_ind_analysis/ind_proj_aligned_fixed/split/'
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
jia.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,jia_traits)
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

## if we compare ERA projections of PC2 with projections for other subtypes is that significant ?
ft <- final[class=='individual' & pc=='PC2' & !trait %in% c('cc','missing','UnA'),]
t.test(ft[trait=='ERA',]$value,mu=final[class=='control' & sample=='ERA_control' & pc=='PC2',]$value)
