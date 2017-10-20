library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)

ml<-list(
  CD = 'bb_CD',
  CEL = 'bb_CEL',
  MS = 'bb_MS',
  RA = 'bb_RA',
  SLE = 'bb_SLE',
  T1D = 'bb_T1D',
  UC = 'bb_UC'
)


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

bb_traits<-unlist(ml)
bb.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,bb_traits)

## leave one out analysis

loo <- lapply(names(ml),function(m){
  message(sprintf("Leaving out %s",m))
  tmp.DT <- basis.DT[basis.DT$trait!=m,]
  ## compute shrinkage using emp
  tmp.shrink.DT<-compute_shrinkage_metrics(tmp.DT)
  M <- create_ds_matrix(tmp.DT,tmp.shrink.DT,'emp')
  ## need to add control where beta is zero
  M <- rbind(M,control=rep(0,ncol(M)))
  pc.emp <- prcomp(M,center=TRUE,scale=FALSE)
  ## project on bb
  bb_M <- create_ds_matrix(bb.DT,tmp.shrink.DT,'emp')
  pred.emp <- predict(pc.emp,newdata=bb_M)
  rbind(pc.emp$x,pred.emp)
})
names(loo)<-names(ml)

g <- function(M){
    M <- cbind(as.data.table(M),trait=rownames(M))
    M$compare<-"none"
    for(i in seq_along(ml)) {
        M[trait %in% c(names(ml)[i], ml[i]), compare:=names(ml)[i]]
    }
    M[trait=="control",compare:="control"]
    M
}

pdf(file='~/tmp/leave_one_out.pdf')
lapply(seq_along(loo),function(i){
  tx<-loo[[i]]
  tx.emp<-g(tx)
  title<-sprintf("Left out %s",names(loo)[i])
  ggplot(tx.emp,aes(x=PC1,y=PC2,color=compare,label=trait,alpha=compare!='none')) +
  geom_point() + theme_bw() + ggtitle(title) + geom_text(show.legend=FALSE) +
  scale_alpha_discrete(guide=FALSE)
})
#ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait,alpha=compare!='none')) +
#geom_point() + theme_bw() + ggtitle('Actual 5913/8829') + geom_text(show.legend=FALSE) +
#scale_alpha_discrete(guide=FALSE)
dev.off()
