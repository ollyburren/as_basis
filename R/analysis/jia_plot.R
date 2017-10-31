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
#m_file<-file.path(support.dir,'as_basis_manifest.tab')
m_file<-file.path(support.dir,'as_basis_manifest_with_jia_cc.tab')
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
bb_traits<-fread(m_file)[grep('jia_',trait),]$trait
#bb_traits<-bb_traits[bb_traits != 'jia_cc']
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
ggplot(emp,aes(x=PC1,y=PC2,color=grepl('jia',trait),label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')

plot <- melt(emp,id.vars=c('trait','compare'))
plot[,label:=ifelse(!grepl('jia',trait),'non.jia',trait)]

## plot scree plots
library(ggplot2)
ggplot(plot,aes(x=variable,y=value,alpha=label!='non.jia',color=label,group=trait)) + geom_point() + geom_line() + theme_bw() + scale_alpha_discrete(range = c(0.1, 1))

#JIA is not adding much I don't think we can probably exclude from the basis.

## worried that MHC is causing the difference so remove chr6 from basis and projection

 no6.basis.DT <- subset(basis.DT, chr !='6')
 no6.shrink.DT<-compute_shrinkage_metrics(no6.basis.DT)
 ## need to add control where beta is zero
no6.basis.mat.emp <- create_ds_matrix(no6.basis.DT,no6.shrink.DT,'emp')
 ## need to add control where beta is zero
  no6.basis.mat.emp<-rbind( no6.basis.mat.emp,control=rep(0,ncol( no6.basis.mat.emp)))
 no6.pc.emp <- prcomp( no6.basis.mat.emp,center=TRUE,scale=FALSE)

no6.emp<-g(no6.pc.emp$x)
ggplot(no6.emp,aes(x=PC1,y=PC2,color=grepl('jia',trait),label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')


no6.bb.DT <- subset(bb.DT, chr !='6')
no6.bb.mat.emp<-create_ds_matrix(no6.bb.DT, no6.shrink.DT,'emp')
no6.pred.emp <- predict( no6.pc.emp,newdata=no6.bb.mat.emp)
no6.emp<-rbind(no6.pc.emp$x,no6.pred.emp)
no6.emp<-g(no6.emp)

ggplot(no6.emp,aes(x=PC1,y=PC2,color=grepl('jia',trait),label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')
