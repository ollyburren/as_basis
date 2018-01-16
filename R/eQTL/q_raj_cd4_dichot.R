library(cupcake)
library(optparse)
DEBUG=FALSE
option_list = list(
        make_option(c("-i", "--idx"), type="character",
                    help="index sequence of probes to compute", metavar="character"),
        make_option(c("-n", "--null"), action="store_true",default=FALSE,
                  help="Flag to switch on simulation under the null", metavar="character"),
        make_option(c("-o", "--out"), type="character",
                  help="Output dir", metavar="character")
        )

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)

  shrink.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/shrink.RDS"
  pc.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/pc.RDS"
  gt.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/gt_matrix.RDS"
  ot.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/outcome_matrix.RDS"
  snps.file="/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/snps.RDS"
  ## this should probably be an argument as not sure we want to overwrite projections as this could get confusing.
  #out.dir <- "/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/projections_null3/"

  probes_seq <- eval(parse(text=args$idx))

  shrink.DT <- readRDS(shrink.file)
  pc.emp <- readRDS(pc.file)
  g.mat <- readRDS(gt.file)
  ot.mat <- readRDS(ot.file)
  snp.ids <- readRDS(snps.file)

if(FALSE){
  s<-seq(1,ncol(ot.mat),by=100)
  range <- data.table(s,c(s[-1]-1,ncol(ot.mat)))[,pid:=paste(s,V2,sep=':')]$pid
  #cmd <- sprintf("Rscript --vanilla /home/ob219/git/as_basis/R/eQTL/q_raj_cd4_dichot.R -i %s",range)
  # this for simulation under the null where we randomly sample to get outcome vector
  cmd <- sprintf("Rscript --vanilla /home/ob219/git/as_basis/R/eQTL/q_raj_cd4_dichot.R -n -i %s -o %s",range,out.dir)
  write(cmd,file="/home/ob219/git/as_basis/sh/raj_cd4_dichot_null3.txt")
}

computeORMat <- function(Y,G,n0=214,n1=212,underNull=FALSE){
  # ploidy is two - two groups are not identical if we use median
  if(is.matrix(Y) & is.matrix(G))
    stop("Both X and Y are matrices in computeORMat")
  ## this code allows us to run things under the null by randomly sampling to
  ## get the out come vector.
  if(underNull){
    #message("Running under null")
    Y<-rep(0L, nrow(G))
    Y[sample.int(length(Y),size=n1/2,replace=FALSE)]<-1
  }
  unx <- colSums((1-Y) * G)
  ex <- colSums(Y * G)
  #(log(n0-unx) + log(ex)) - (log(unx) + log(n1-ex)) ##==
  log(n0-unx) + log(ex) - log(unx) - log(n1-ex)
}

ret <- lapply(probes_seq,function(pidx){
  pout <- ot.mat[,pidx]
  computeORMat(pout,g.mat,underNull=args$null)
}) %>% do.call("rbind",.)

## as input expect OR

ret <- exp(ret) %>% t %>% data.table
ret[,snp:=snp.ids]
setnames(ret,c(as.character(colnames(ot.mat)[probes_seq]),'pid'))
DT <- melt(ret,id.vars='pid')
setnames(DT,c('pid','sample','or'))
## some OR are inf assume because missing allele ?
DT[is.infinite(DT$or) | DT$or==0 ,or:=1]

setcolorder(DT,c('sample','pid','or'))
setkey(DT,pid)
setnames(DT,'sample','trait')
mat.emp <- create_ds_matrix(DT,shrink.DT,'emp')
proj<-predict(pc.emp,newdata=mat.emp)
saveRDS(proj,file.path(args$out,paste(args$idx,'RDS',sep='.')))

if(FALSE){
  library(ggplot2)
  all.proj.null <- lapply(list.files(path='/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/projections_null/',pattern="*.RDS",full.names=TRUE),readRDS)
  all.proj.null <- do.call('rbind',all.proj.null)
  all.DT.null <- data.table(all.proj.null)[,c('probe','type'):=list(rownames(all.proj.null),'null')]

  all.proj <- lapply(list.files(path='/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/projections/',pattern="*.RDS",full.names=TRUE),readRDS)
  all.proj <- do.call('rbind',all.proj)
  all.DT <- data.table(all.proj)[,c('probe','type'):=list(rownames(all.proj),'real')]

  all.DT <- rbind(all.DT.null,all.DT)

  #all.DT <- melt(all.proj,id.vars=probe)
  library(cowplot)
  ggplot(all.DT,aes(x=PC1,y=PC2)) + geom_point(alpha=0.2,size=0.3) +
  geom_point(data=data.table(pc.emp$x)[,label:=rownames(pc.emp$x)],color='red',size=1) + geom_text(data=data.table(pc.emp$x)[,label:=rownames(pc.emp$x)],aes(label=label,color='red')) + facet_wrap(~type)


  ggplot(all.DT,aes(x=PC2,y=PC4)) + geom_point(alpha=0.2,size=0.3) +
  geom_point(data=data.table(pc.emp$x)[,label:=rownames(pc.emp$x)],aes(x=PC2,y=PC4),color='red',size=1) + geom_text(data=data.table(pc.emp$x)[,label:=rownames(pc.emp$x)],aes(label=label,color='red'))

  ##  check top genes from individual analysis
  tgenes <- readRDS("~/tmp/top_genes.RDS")

  pc.cmp <- 'PC4'
  ## get PC's for top probes
  plots <- lapply(paste0('PC',c(2,4,7)),function(pc.cmp){
    topy <- all.DT[paste('P',probe,sep='_') %in% tgenes[PC==pc.cmp,]$lu]
    ggplot(all.DT,aes_q(as.name(pc.cmp))) + geom_histogram() + geom_vline(data=topy,aes_q(xintercept=as.name(pc.cmp)),color='firebrick') + facet_wrap(~type)
  })
  pg <- plot_grid(plots[[1]],plots[[2]],plots[[3]],labels=paste0('PC',c(2,4,7)))


  ## compute p.vals for each of our top 10 genes by the number of genes with a more extreme p.val

res.p <- lapply(paste0('PC',1:10),function(pc.cmp){
  topy <- all.DT[paste('P',probe,sep='_') %in% tgenes[PC==pc.cmp,]$lu & type=='real'][[pc.cmp]]
  prbs <- all.DT[paste('P',probe,sep='_') %in% tgenes[PC==pc.cmp,]$lu & type=='real'][['probe']]
  ctrl <-  pc.emp$x["control",pc.cmp]
  ## distances
  d <- abs(ctrl-all.proj.null[,pc.cmp])
  tg <- abs(ctrl - topy)
  p.vals <- sapply(tg,function(x){
    sum(x<d)/length(d)
  })
  p.vals[p.vals==0] <- 1/length(d)
  names(p.vals) <- paste('P',prbs,sep='_')
  data.table(ind.mlp=-log10(tgenes[PC==pc.cmp,]$p.coeff),agg.mlp=-log10(p.vals[tgenes[PC==pc.cmp,]$lu]),pc=pc.cmp)
})

all.cmp<-rbindlist(res.p)[,pc:=factor(pc,levels=paste0('PC',1:10))]
ggplot(all.cmp,aes(x=ind.mlp,y=agg.mlp)) + geom_point() + facet_wrap(~pc) + xlab("Individual -log10(P)") + ylab("Aggregate -log10(P)") + geom_abline(intercept=0,slope=1,lty=2,color="steelblue")

agenes <- readRDS("/home/ob219/scratch/as_basis/raj_cd4_ind_analysis/regression_analysis/all_genes.RDS")

es.p <- lapply(paste0('PC',1:10),function(pc.cmp){
  topy <- all.DT[paste('P',probe,sep='_') %in% agenes[PC==pc.cmp,]$lu & type=='real'][[pc.cmp]]
  prbs <- all.DT[paste('P',probe,sep='_') %in% agenes[PC==pc.cmp,]$lu & type=='real'][['probe']]
  ctrl <-  pc.emp$x["control",pc.cmp]
  ## distances
  d <- abs(ctrl-all.proj.null[,pc.cmp])
  tg <- abs(ctrl - topy)
  p.vals <- sapply(tg,function(x){
    sum(x<d)/length(d)
  })
  p.vals[p.vals==0] <- 1/length(d)
  names(p.vals) <- paste('P',prbs,sep='_')
  data.table(ind.mlp=-log10(agenes[PC==pc.cmp,]$p.coeff),agg.mlp=-log10(p.vals[agenes[PC==pc.cmp,]$lu]),pc=pc.cmp)
})

saveRDS(es.p,file="/scratch/ob219/tmp/eqtl.pvals.RDS")

all.acmp<-rbindlist(es.p)[,pc:=factor(pc,levels=paste0('PC',1:10))]
ggplot(all.acmp,aes(x=ind.mlp,y=agg.mlp)) + geom_hex() + facet_wrap(~pc) + xlab("Individual -log10(P)") + ylab("Aggregate -log10(P)") + geom_abline(intercept=0,slope=1,lty=2,color="steelblue")

es.p

}
