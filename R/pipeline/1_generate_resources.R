  source('~/git/as_basis/R/pipeline/as_basis_functions.R')

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
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat)))
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
  ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw()

  ## what happens if we use estimate instead ?
  basis.mat.est <- create_ds_matrix(basis.DT,shrink.DT,'est')
  ## need to add control where beta is zero
  basis.mat.est<-rbind(basis.mat.est,control=rep(0,ncol(basis.mat)))
  pc.est <- prcomp(basis.mat.est,center=TRUE,scale=FALSE)
  bb.mat.est<-create_ds_matrix(bb.DT,shrink.DT,'est')
  pred.emp <- predict(pc.est,newdata=bb.mat.est)
  est<-rbind(pc.emp$x,pred.emp)

  est<-g(est)
  ggplot(est,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw()



  if(FALSE){
    ## plot to check things
    library(ggplot2)
    ## take a random sample of snps
    keep<-sample(gh_maf.DT$pid,10000)
    toplot<-data.table::melt(unique(basis.DT[,.(pid,maf,mean_maf_se,est_maf_se)],by='pid'),id.vars=c('pid','maf'),measure.vars=c('mean_maf_se','est_maf_se'))
    sdata<-toplot[toplot$pid %in% keep,]
    ggplot(sdata,aes(x=maf,y=value,color=variable)) + geom_point()
  }
