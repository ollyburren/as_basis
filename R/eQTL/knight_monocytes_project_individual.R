library(cupcake)
support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support_cd14_eQTL.tab')
m_file<-file.path(support.dir,'as_manifest_december.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'

basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)


input.dir <- '/scratch/ob219/as_basis/raj_cd14_ind_analysis/eqtl_ind_proj_2_0.01_2500_1e+06'
in.files <- list.files(path=input.dir,pattern='*.RDS',full.names=TRUE)

library(parallel)
all.preds <- mclapply(in.files,function(f){
  message(f)
  obj <- readRDS(f)
  snps <- obj$snps
  or <- data.table(exp(obj$proj.lor))
  or[,snp:=snps$pid]
  samples <- obj$samples
  setnames(or,c(as.character(samples$pedigree),'pid'))
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
all.proj <- dcast(all.proj,sample~pc)
library(ggplot2)
library(cowplot)
vexp <- summary(pc.emp)[['importance']][2,]
PC1.var<-signif(vexp["PC1"]*100,digits=3)
PC2.var<-signif(vexp["PC2"]*100,digits=3)
ppr<-ggplot(all.proj,aes(x=PC1,y=PC2)) + geom_point(size=3) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none")

## next we need to regress PC's vs gene expression

## load in gene expression values

(load("/scratch/wallace/twas/raj-cd14-expression.RData"))
expr.DT<-data.table(t(expr))
setnames(expr.DT,paste('P',rownames(expr),sep='_'))
expr.DT[,sample:=colnames(expr)]
## load the transform matrix
translate <- readRDS('/scratch/wallace/twas/model_output/trans_raj-cd14.rds')
expr.DT[,sample_id:=translate[sample]]
setkey(expr.DT,sample_id)
setkey(all.proj,sample)
for.reg <- all.proj[expr.DT]

ggplot(for.reg,aes(x=PC1,y=P_7896908)) + geom_point() + geom_smooth(method='lm')
saveRDS(for.reg,file="/home/ob219/scratch/as_basis/raj_cd14_ind_analysis/regression_analysis/cd4.RDS")
for.reg <- readRDS("/home/ob219/scratch/as_basis/raj_cd14_ind_analysis/regression_analysis/cd4.RDS")


## first we attempt to fit a set of linear models where PC is the explantory model
library(parallel)
all.lms <- mclapply(names(for.reg)[grep("PC",names(for.reg))],function(pc){
  message(pc)
  lapply(grep('^P\\_',names(for.reg)),function(i){
    probe <- names(for.reg)[i]
    dat <- for.reg[,c(pc,probe),with=FALSE]
    setnames(dat,c('PC','expression'))
    lm(expression~PC,dat)
  })
},mc.cores=8)

names(all.lms)<-names(for.reg)[grep("PC",names(for.reg))]
## for each model we can obtain t-statistics for the beta coefficients as p-values.
## we can the check for inflation by plotting as a qqplot trellis plot.

## should create a library for this.
get_qq_dt <- function(x,l=0.99,minx=TRUE){
  n <- length(x)
  ## expected
  q <- -log10((n:1)/(n+1))
  ## observed
  x <- sort(x)
  ## l is a parameter that allows us to only plot a subset of the data
  n1 <- round(l*n)
  if(minx)
    return(data.table(expected=c(0,q[n1:n]),observed=c(0,x[n1:n])))
  return(data.table(expected=q,observed=x))
}

pc.p <- lapply(seq_along(names(all.lms)),function(i){
  pc<-names(all.lms)[i]
  all.p <- sapply(all.lms[[i]],function(x){
    summary(x)$coefficient["PC",4]
  })
  DT <- data.table(PC=pc,probe=names(for.reg)[grep('^P\\_',names(for.reg))],p.coeff=all.p)
})
## for qq.plot

for.qq <- rbindlist(lapply(pc.p,function(x) cbind(unique(x$PC),get_qq_dt(-log10(x$p.coeff)))))
for.qq$V1<-factor(for.qq$V1,levels=paste0('PC',1:10))
## perhaps hilight certain PC's
hl.pc<-paste0('PC',c(20))
for.qq[,hilight:=V1 %in% hl.pc]

library(cowplot)
gpl <- ggplot(for.qq,aes(x=expected,y=observed,color=hilight)) +
geom_abline(intercept=0,slope=1,color='dodgerblue',lty=2) +
geom_point(size=0.5) + facet_wrap(~V1,ncol=3,nrow=4) +
ylab("Observed -log10(P)") + xlab("Expected -log10(P)") +
scale_color_manual(guide=FALSE,values=c('TRUE'='firebrick','FALSE'='black'))
save_plot("~/tmp/raj_cd14_qq.pdf",gpl,base_height=5)

## Chris wants top genes for PC's from qq plots looks as if a threshold of 10^-3

## annotate probes with genes

raj.anno <- data.table(readRDS("/scratch/wallace/twas/gene_details.rds"))[,lu:=paste('P',affy_hugene_1_0_st_v1,sep='_')]
setkey(raj.anno,lu)

anno<-lapply(pc.p,function(x){
  setkey(x,probe)
  tmp<-raj.anno[x]
  dup.probes <- tmp[duplicated(tmp$lu),]$lu
  tmp<-tmp[!lu %in% dup.probes,]
  head(tmp[order(-log10(p.coeff),decreasing=TRUE),],n=10)
})

saveRDS(rbindlist(anno),"/home/ob219/scratch/as_basis/raj_cd14_ind_analysis/regression_analysis/top_genes.RDS")

## do all genes as well for comparison

anno.all<-lapply(pc.p,function(x){
  setkey(x,probe)
  tmp<-raj.anno[x]
  dup.probes <- tmp[duplicated(tmp$lu),]$lu
  tmp<-tmp[!lu %in% dup.probes,]
  tmp[order(-log10(p.coeff),decreasing=TRUE),]
})

saveRDS(rbindlist(anno.all),"/home/ob219/scratch/as_basis/raj_cd14_ind_analysis/regression_analysis/all_genes.RDS")
