## This script simulates what happens under the null i.e when the scaled beta's i.e. gamma hats are equal to one another.

## assumes that cupcake has been loaded with (devtools)
library(cupcake)
library(optparse)

TEST <- FALSE # set to true when debugging
DEFAULT.SUPPORT.DIR <- '/home/ob219/rds/hpc-work/as_basis/support_tab'
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_feb/'
DEFAULT.N.SIMS <- 50
DEFAULT.GWAS.DATA.DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned/'
#DEFAULT.MANIFEST.FILE <- 'as_basis_manifest_with_jia_cc.tab'
DEFAULT.MANIFEST.FILE <- '/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv'
DEFAULT.REF.AF.FILE <- '/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab';
## does not matter what this is as long as has correct number of SNPs !
DEFAULT.TRAIT <- 'ret_p'
#DEFAULT_CACHE_BASIS <- '/home/ob219/rds/hpc-work/as_basis/cache/as_basis_snp_support_feb_2018_w_ms.RDS'
DEFAULT_CACHE_BASIS <- '/home/ob219/rds/hpc-work/as_basis/cache/as_basis_snp_support_feb_2018_w_ms_ws.RDS'
#SHRINKAGE_METHOD <- 'emp'
SHRINKAGE_METHOD <- 'ws_emp'

## this routine takes a set of vcf files from a reference panel say 1KG and extracts the variants that are used in our basis.
#In this form it won't work (needs basis info for filtering but acts as a reference)

filterVCF <- function(){
  study1.DT<-get_gwas_data(DEFAULT.MANIFEST.FILE,DEFAULT.REF.AF.FILE,DEFAULT.GWAS.DATA.DIR,FALSE,DEFAULT.TRAIT )
  study1.DT[,c('chr','position'):=tstrsplit(pid,':')]
  bcft <- '~/bin/bcftools-1.4/bcftools'
  vcf.dir<-'/home/ob219/rds/hpc-work/DATA/1kgenome/VCF/EUR/by.chr.phase1/'
  out.dir <- DEFAULT.SNPSTATS.DIR
  s.DT <- split(study1.DT[,.(chr,position)],study1.DT$chr)
  for(chr in names(s.DT)){
      message(sprintf("Processing %s",chr))
      tmp<-s.DT[[chr]][order(as.numeric(chr),as.numeric(position)),]
      tmp.file <- tempfile(pattern = "vcf_tmp", tmpdir = tempdir())
      write.table(tmp,file=tmp.file,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      #vcf.file <- file.path(vcf.dir,sprintf("ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz",chr))
      vcf.file <- file.path(vcf.dir,sprintf("chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz",chr))
      obj <- vcf2snpmatrix(vcf.file,bcft,tmp.file,TRUE)
      unlink(tmp.file)
      out.file <- file.path(out.dir,sprintf("%s_1kg.RData",chr))
      save(obj,file=out.file)
    }
}

## to prevent unneccsary computation we can compute the basis ahead of time and just read that in

cacheBasis <- function(){
  basis.DT<-get_gwas_data(DEFAULT.MANIFEST.FILE,DEFAULT.REF.AF.FILE,DEFAULT.GWAS.DATA.DIR,FALSE)
  shrink.DT<-compute_shrinkage_metrics(basis.DT)
  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  saveRDS(list(basis=pc.emp,shrink=shrink.DT),file=DEFAULT_CACHE_BASIS)
}

## perhaps add options for different manifests etc.
option_list = list(
        make_option(c("-n","--total"), type="numeric",default=DEFAULT.SUPPORT.DIR,
                    help="Location of support files", metavar="character"),
        make_option(c("-c", "--ncases"), type="numeric", default=DEFAULT.GWAS.DATA.DIR,
                    help="location of OR aligned GWAS source files", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                    help="output directory", metavar="character")
        )


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  print(args)
}else{
  message("IN TESTING MODE ======>!")
  trait <- 'ret_p'
  args <- list(
    total = 10000,
    ncases=5000,
    out_dir = '~/tmp/test_distance/'
  )
  print(args)
}


## simulate
cache.obj <- readRDS(DEFAULT_CACHE_BASIS)
study1.DT<-get_gwas_data(DEFAULT.MANIFEST.FILE,DEFAULT.REF.AF.FILE,DEFAULT.GWAS.DATA.DIR,FALSE,DEFAULT.TRAIT )
study1.DT <- study1.DT[cache.obj$shrink[,.(pid,emp_maf_se)]]
## simulate under the null
## need tp compute empirical SE for study
n1 <- args$ncases
n <- args$total
# note that emp_maf_se

## compute the standard error associatied with our


est_SE<-function(f,n.cases,total){
  n.ctrls<-total-n.cases
  a<-1/(2*f*n.ctrls)
  b<-1/(2*(1-f)*n.ctrls)
  c<-1/(2*f*n.cases)
  d<-1/(2*(1-f)*n.cases)
  sqrt(rowSums(cbind(a,b,c,d)))
  #return(sqrt(sum(ca,b,c,d))))
}

if(FALSE){
  ## work out scale factor so we don't need to compute est_SE for each sample size
  #number of SNPs
  n<-10
  mafs<-runif(n)
  tot <- 2000
  n.cases1<-1000
  n.cases2<-500
  n.controls1<-tot-n.cases1
  n.controls2<-tot-n.cases2
  se1 <- est_SE(mafs,n.cases1,total)
  se2 <- est_SE(mafs,n.cases2,total)
  cfact <- unique(se1/se2)[1]
  prop1 <- n.cases1/total
  prop2 <- n.cases2/total
  dfact <- sqrt(1/((prop * (1 - prop)) / (prop2 * (1-prop2))))
  identical(dfact,cfact)
  se1cp <- tcrossprod(se1)
  se2cp <- tcrossprod(se2)
  # make up some weights
  w<-runif(10)
  ## make up some sigma for n
  A <- matrix(runif(n^2)*2-1, ncol=n)
  Sigma1 <- tcrossprod(A) * tcrossprod(se1)
  Sigma2 <- tcrossprod(A) * tcrossprod(se2)
  c1<-sum(tcrossprod(w) * Sigma1)
  c2<-sum(tcrossprod(w) * Sigma2)
  multi <-  (prop2 * (1-prop2))/(prop * (1-prop))
  c2prime <- sum(tcrossprod(w) * Sigma1) * multi
}


study1.DT$n <- NULL
study1.DT$n1 <- NULL

#study1.DT[,emp_se := 1/sqrt(2) * sqrt(1/n1 + (1/sqrt(n-n1))) * emp_maf_se]
#n<-8000
#n1<-4000
study1.DT[,emp_se_old := 1/sqrt(2) * sqrt(1/n) * emp_maf_se]
#study1.DT[,est_se:=est_SE(maf,n1,n)]
study1.DT[,emp_se:=est_SE(maf,n1,n)]
#study1.DT[,emp_se_old := 1/sqrt(2) * 1/sqrt(n) * emp_maf_se]
study1.DT[,or:=1]
## need chr and perhaps position
study1.DT[,c('chr','position'):=tstrsplit(pid,':')]
message(sprintf("Simulate %d studies",args$n_sims))
s1.sim <- simulate_study(study1.DT,DEFAULT.SNPSTATS.DIR,shrink_beta=FALSE,n_sims=DEFAULT.N.SIMS,quiet=TRUE)
## compute the projection on 10 at a time to prevent huge matrix generation eating all memory
idx <- split(1:nrow(s1.sim),s1.sim$trait)
pred.emp<-do.call('rbind',lapply(split(idx, ceiling(seq_along(idx)/10)),function(i){
  tmp.idx <- unlist(i)
  tmp.sim <- s1.sim[tmp.idx,]
  setkey(tmp.sim,pid)
  sim.mat.emp<-create_ds_matrix(tmp.sim,cache.obj$shrink,SHRINKAGE_METHOD)
  predict(cache.obj$basis,newdata=sim.mat.emp)
}))
rand.string <- tolower(do.call(paste0, replicate(6, sample(LETTERS, 1, TRUE), FALSE)))
of <- sprintf("%d_%d_%s.RDS",args$total,args$ncases,rand.string)
ofile <- file.path(args$out_dir,of)
saveRDS(pred.emp,file=ofile)
message(sprintf("Saved projections in %s",ofile))


if(FALSE){
  library(data.table)
  library(magrittr)
  library(cowplot)
  OUTDIR <- '/home/ob219/rds/hpc-work/as_basis/callibration/filter_feb_2018_w_ms_ws_fixed'
  cmd_template<-"Rscript /home/ob219/git/as_basis/R/analysis/null_by_sample_size_special.R -n %d -c %d -o %s"
  all.cmds <- lapply(c(350000,15000,5000),function(n){
    sapply(c(250,500,1000,5000,10000,50000),function(n1){
      if(n1<n)
        rep(sprintf(cmd_template,n,n1,OUTDIR),40)
    })
  }) %>% unlist
  write(all.cmds,"~ob219/git/as_basis/sh/run_cal.txt")

  all.cmds <- lapply(c(1000,2000,3000),function(n){
    sapply(c(250,500,1000,5000,10000,50000),function(n1){
      if(n1<n)
        rep(sprintf(cmd_template,n,n1,OUTDIR),40)
    })
  }) %>% unlist
  write(all.cmds,"~ob219/git/as_basis/sh/run_cal2.txt")

  ## compile callibrations

  af <- list.files(path=OUTDIR,pattern="*.RDS",full.names=TRUE)
  byexp<-split(af,gsub("\\_[a-z]+.RDS$","",basename(af)))

  ## check for normality
  all.res <- lapply(byexp[['5000_250']],readRDS) %>% do.call('rbind',.)
  qqnorm(all.res[,"PC1"])
  qqline(all.res[,"PC1"])
  ## do wilcox testing
  test<-all.res[,"PC1"]
  ks.test((test-mean(test))/sd(test),"pnorm")
  qqnorm((test-mean(test))/sd(test))

  ## build callibration curve for one sample size to start to see if there is a trend

  all.res <- lapply(names(byexp),function(x){
    all.res <- lapply(byexp[[x]],readRDS) %>% do.call('rbind',.)
    var <- apply(all.res,2,var)
    tmp<-strsplit(x,"\\_") %>% unlist
    data.table(total=tmp[1],ncases=tmp[2],variance=var,pc=names(var))
  }) %>% rbindlist


  ## if our estimate of the variance is drawn
  CQF.lt <- qchisq(0.025, 500-1, lower.tail=FALSE)
  CQF.ut <- qchisq(0.025, 500-1, lower.tail=TRUE)

  ## an approximate way to do this is variance * sqrt(2/(n-1)) ## I should ask Chris how this works



  all.res[,c('ci.lower','ci.upper'):=list((variance * (500-1)/CQF.lt),((variance * (500-1))/CQF.ut)) ]
  library(scales)
  all.res[,ci:=variance * sqrt(2/(500-1))]
  all.res[,total:=factor(total,levels=sort(unique(as.numeric(total))))]
  pd <- position_dodge(width=0.3)
  # ppf <- ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
  # geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + xlab("Cases") + scale_color_discrete("Sample Size") + ylab("Variance (PC1 Loadings)") +
  # scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + theme(legend.position =c(0.65,0.85)) +
  # background_grid(major = "xy", minor = "none")
  #
  #
  # all.res[,total:=as.numeric(as.character(total))]
  # all.res[,ncases:=as.numeric(as.character(ncases))]
  # library(magrittr)
  # foo<-lapply(paste0('PC',1:10),function(PC){
  #     mod<-summary(lm(data=all.res[pc==PC,],log(variance)~log(total) + log(ncases)))$coefficients
  #     data.table(pc=PC,name=rownames(mod),mod)
  #   }) %>% rbindlist





  ## used AIC to get the model

  # model that we can use is

  calli<-lm(data=all.res[pc!="PC11",],log(variance)~log(total)+log(ncases) + pc)
  all.res[,prop1:=ncases/total]
  all.res[,prop2:=(1-prop1) * prop1]
  saveRDS(calli,"/home/ob219/rds/hpc-work/as_basis/callibration/callibration_model_fixed.RDS")
  ## do the

  ## attempt things analytically
  ## simulate
#
#   compute_proj_var <- function(DT,w.DT,shrink.DT,ref_gt_dir,method='ws_emp',quiet=TRUE){
#   s.DT <- split(DT,DT$chr)
#   setkey(w.DT,pid)
#   all.chr <- lapply(names(s.DT),function(chr){
#     if(!quiet)
#       message(sprintf("Processing %s",chr))
#     ss.file<-file.path(ref_gt_dir,sprintf("%s_1kg.RData",chr))
#     sm<-get(load(ss.file))
#     ## there are sometimes duplicates that we need to remove
#     dup.idx<-which(duplicated(obj$info$pid))
#     if(length(dup.idx)>0){
#       if(!quiet)
#         message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
#       sm$info<-sm$info[-dup.idx,]
#       sm$sm <- sm$sm[,-dup.idx]
#     }
#     sm$info$order<-1:nrow(sm$info)
#     # by ld block
#     by.ld <- split(s.DT[[chr]],s.DT[[chr]]$ld.block)
#     chr.var <- lapply(names(by.ld),function(block){
#       if(!quiet)
#         message(sprintf("Processing %s",block))
#       dat <- by.ld[[block]]
#       w <- w.DT[pid %in% dat$pid,]
#       vmethod = sprintf("%s_shrinkage",method)
#       s <- shrink.DT[pid %in% dat$pid,c('pid',vmethod),with=FALSE]
#       setkey(w,pid)
#       setkey(dat,pid)
#       dat <- dat[w]
#       setkey(s,pid)
#       dat <- dat[s]
#       info <-sm$info[pid %in% dat$pid ,.(pid,order)]
#       setkey(info,pid)
#       dat <- dat[info][order(order),]
#       ## need to add shrinkage if we attempt to compute under the alternative
#       Sigma <- with(dat,cov_beta(sm$sm[,order],emp_se))
#       pc.cols <- which(grepl("^PC[0-9]+$",names(dat)))
#       # multiply basis weights by the shrinkge metric
#       w.mat <- as.matrix(dat[,pc.cols,with=FALSE]) * dat[[vmethod]]
#       #pc^{T} %*% pc * Sigma
#       apply(w.mat,2,function(pc) sum(tcrossprod(pc) * Sigma))
#     })
#     colSums(do.call('rbind',chr.var))
#   })
#   colSums(do.call('rbind',all.chr))
# }


  library(cupcake)
  cache.obj <- readRDS(DEFAULT_CACHE_BASIS)
  study1.DT<-get_gwas_data(DEFAULT.MANIFEST.FILE,DEFAULT.REF.AF.FILE,DEFAULT.GWAS.DATA.DIR,FALSE,DEFAULT.TRAIT )
  study1.DT <- study1.DT[cache.obj$shrink[,.(pid,emp_maf_se)]]
  w.DT <- data.table(pid=rownames(cache.obj$basis$rotation),cache.obj$basis$rotation)
  shrink.DT <- cache.obj$shrink

  confs <- unique(all.res[,.(total,ncases)])
  study1.DT[,c('chr','position'):=tstrsplit(pid,':')]

  ## this takes a while to run over the 22 configurations.

  vars <- lapply(1:nrow(confs),function(i){
    message(i)
    study1.DT[,n1:=confs$ncases[i]]
    study1.DT[,n:=confs$total[i]]
    #study1.DT[,c('n','n1'):=list(n,n1)]
    study1.DT[,emp_se:=est_SE(maf,n1,n)]
    compute_proj_var(study1.DT,w.DT,cache.obj$shrink,DEFAULT.SNPSTATS.DIR,'ws_emp',quiet=TRUE)
  })

  saveRDS(vars,file="~/tmp/vars.RDS")

  all.vars <- do.call('rbind',vars)
  all.vars <- cbind(confs,all.vars)
  all.vars <- melt(all.vars,id.vars=c('total','ncases'),measured.vars=paste0('PC',1:11))
  setnames(all.vars,c('total','ncases','pc','variance'))
  all.vars[,total:=factor(total,levels=levels(all.res$total))]

  ppf <- ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
  geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + xlab("Cases") + scale_color_discrete("Sample Size") +
  ylab("Variance (PC1 Loadings)") + scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position =c(0.65,0.8)) + background_grid(major = "xy", minor = "none") +
  geom_point(data=all.vars[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,color=total,fill=total),position=pd,pch=17,size=3,alpha=0.4,inherit.aes=FALSE) +
  scale_fill_discrete(guide=FALSE)

  save_plot("~/tmp/var_simulation_null.pdf",ppf,base_width=7)

  ## what about scaling ?

  all.vars[,prop:=ncases/as.numeric(as.character(total))]
  all.vars[,factor:=prop * (1-prop)]

  ## derive the variance for all others using 350000
  r.DT <- subset(all.vars, total==15000 & ncases==5000)
  rf.factor <- r.DT$factor
  rf.variance <- r.DT$variance
  all.vars[,dvar:=sqrt(rf.factor/factor) * variance]

  ppf <- ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
  geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + xlab("Cases") + scale_color_discrete("Sample Size") +
  ylab("Variance (PC1 Loadings)") + scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position =c(0.65,0.85)) + background_grid(major = "xy", minor = "none") +
  geom_point(data=all.vars[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=dvar,color=total,fill=total),position=pd,pch=17,size=5,alpha=0.4,inherit.aes=FALSE) +
  scale_fill_discrete(guide=FALSE)

  ## also need to add boot strap points ?

}
