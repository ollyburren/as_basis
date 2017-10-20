library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)

  ## import GWAS data for basis

  ## support files
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
  ## need to add control where beta is zero
  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

  ## project on biobank to see if we can recreate Chris' figure.

  bb_traits<-fread(m_file)[grep('bb_',trait),]$trait
  bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,bb_traits)
  bb.mat.emp<-create_ds_matrix(bb.DT,shrink.DT,'emp')
  pred.emp <- predict(pc.emp,newdata=bb.mat.emp)

  emp<-rbind(pc.emp$x,pred.emp)

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
  ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')

  ## what happens if we use estimate instead ?
  basis.mat.est <- create_ds_matrix(basis.DT,shrink.DT,'est')
  ## need to add control where beta is zero
  basis.mat.est<-rbind(basis.mat.est,control=rep(0,ncol(basis.mat.est)))
  pc.est <- prcomp(basis.mat.est,center=TRUE,scale=FALSE)
  bb.mat.est<-create_ds_matrix(bb.DT,shrink.DT,'est')
  pred.emp <- predict(pc.est,newdata=bb.mat.est)
  est<-rbind(pc.emp$x,pred.emp)
  est<-g(est)
  ggplot(est,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Estimate MAF SE shrinkage')


## This is code to compare the MAF shrinkage.
  if(FALSE){
    ## plot to check things
    library(ggplot2)
    ## take a random sample of snps
    #keep<-sample(shrink.DT$pid,10000)
    toplot <- shrink.DT[unique(basis.DT[basis.DT$pid %in% sample(unique(shrink.DT$pid),10000),.(maf),by=pid])]
    toplot[,ratio:=emp_maf_se/est_maf_se]
    ggplot(toplot,aes(x=maf,y=ratio)) + geom_point() + theme_bw()
    toplot<-data.table::melt(toplot[,.(pid,maf,emp_maf_se,est_maf_se)],id.vars=c('pid','maf'),measure.vars=c('emp_maf_se','est_maf_se'))
    ggplot(toplot,aes(x=maf,y=value,color=variable)) + geom_point() + theme_bw()
  }

  ## compute eucledian distance pairwise
pairwise_euc_matrix(est,summary(pc.est)$importance[2,]) %>% dist %>% hclust %>% plot
pairwise_euc_matrix(est,summary(pc.est)$importance[2,],1) %>% dist %>% hclust %>% plot

## next we want to titrate in different case control sizes and

sizes<-c(50,100,200,800,1600,3000,4500,6000)
## for each one need to recompute shrinkages anyways
sim.DT<-basis.DT[basis.DT$trait!='T1D',]
sims<-lapply(sizes,function(s){
  message(sprintf("Processing %d",s))
  tmp.DT<-shrink.DT[basis.DT[basis.DT$trait=='T1D',]]
  ## we assume that OR estimate is constant the thing that changes is it's standard error
  ## we already computed an estimate of standard error due to MAF across all basis use this
  ## for time being but perhaps get for specific trait to compare ?
  ## se(theta) = 1/sqrt(2) * 1/sqrt(n) * sd.maf
  se.lor <- 1/sqrt(2) * 1/sqrt(2*s) * tmp.DT$emp_maf_se
  ## compute epsilon which is a noise parameter to modulate
  ## the OR generated.
  lor <- log(tmp.DT$or) + rnorm(1,0,se.lor)
  tmp.DT[,c('or','p.val','n','n1') := list(exp(lor),2* pnorm(abs(lor/se.lor)),2*s,s)]
  #tmp.DT[,p.val := 2* pnorm(abs(log(or)/se.lor), lower.tail = FALSE)]
  nb.DT <- rbind(sim.DT,tmp.DT[,names(sim.DT),with = FALSE])
  setkey(nb.DT,pid)
  nb.shrink.DT<-compute_shrinkage_metrics(nb.DT)
  ## need to add control where beta is zero
  nb.mat.emp <- create_ds_matrix(nb.DT,nb.shrink.DT,'emp')
  ## need to add control where beta is zero
  nb.mat.emp<-rbind(nb.mat.emp,control=rep(0,ncol(nb.mat.emp)))
  pc.nb <- prcomp(nb.mat.emp,center=TRUE,scale=FALSE)
  nb.bb.mat.emp<-create_ds_matrix(bb.DT,nb.shrink.DT,'emp')
  pred.emp <- predict(pc.nb,newdata=nb.bb.mat.emp)
  list(pc=pc.nb,shrink=copy(nb.shrink.DT),bb.pred=pred.emp)
})

## quick look using PC1 and PC2

pdf(file='~/tmp/t1d_titration_epsilon.pdf')
lapply(seq_along(sims),function(i){
  sa<-sims[[i]]
  tx<-rbind(sa$pc$x,sa$bb.pred)
  tx.emp<-g(tx)
  title<-sprintf("T1D Titration %d cases",sizes[i])
  ggplot(tx.emp,aes(x=PC1,y=PC2,color=compare,label=trait,alpha=compare!='none')) +
  geom_point() + theme_bw() + ggtitle(title) + geom_text(show.legend=FALSE) +
  scale_alpha_discrete(guide=FALSE)
})
ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait,alpha=compare!='none')) +
geom_point() + theme_bw() + ggtitle('Actual 5913/8829') + geom_text(show.legend=FALSE) +
scale_alpha_discrete(guide=FALSE)
dev.off()
stop()



## create basis without T1D and project back on

not1d.DT<-basis.DT[basis.DT$trait!='T1D',]
not1d.shrink.DT<-compute_shrinkage_metrics(not1d.DT)
## need to add control where beta is zero
not1d.mat.emp <- create_ds_matrix(not1d.DT,not1d.shrink.DT,'emp')
## need to add control where beta is zero
not1d.mat.emp<-rbind(not1d.mat.emp,control=rep(0,ncol(not1d.mat.emp)))
pc.not1d <- prcomp(not1d.mat.emp,center=TRUE,scale=FALSE)
t1d.bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,c(bb_traits,'T1D'))
not1d.bb.mat.emp<-create_ds_matrix(t1d.bb.DT,not1d.shrink.DT,'emp')
no.t1d.pred.emp <- predict(pc.not1d,newdata=not1d.bb.mat.emp)
tx<-rbind(pc.not1d$x,no.t1d.pred.emp)
tx.emp<-g(tx)
title<-"NO T1D"
ggplot(tx.emp,aes(x=PC1,y=PC2,color=compare,label=trait,alpha=compare!='none')) +
geom_point() + theme_bw() + ggtitle(title) + geom_text(show.legend=FALSE) +
scale_alpha_discrete(guide=FALSE)
