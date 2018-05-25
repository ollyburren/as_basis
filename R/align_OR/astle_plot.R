library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)
## import GWAS data for basis
## support files
support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support.tab')
#ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-file.path(support.dir,'as_manifest_december.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
## project on biobank to see if we can recreate Chris' figure.
bb_traits<-fread(m_file)[grep('bb_',trait),]$trait
bb_traits<-c(bb_traits,'breast_cancer')
bb.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,bb_traits)
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
  UC = 'bb_UC',
  PBC = 'PBC',
  PSC = 'PSC'
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
emp<-g(emp)[,is.biobank:=FALSE]
emp[grep("^bb",trait),is.biobank:=TRUE]
emp[,is.imb.ai:=compare %in% names(ml)]
emp[,glabel:=trait]
emp[is.imb.ai==FALSE,glabel:='']
emp[trait=='control',glabel:='CONTROL']
emp[trait=='bb_colitis',glabel:='bb_colitis_not_crohns_or_ulcerative_colitis']
#ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')

PC1.var<-signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3)
PC2.var<-signif(summary(pc.emp)[['importance']][2,]["PC2"]*100,digits=3)
## nice plot
library("cowplot")
ppl<-ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=glabel,alpha=glabel!='')) + geom_point(size=3) + geom_text(hjust = 0, nudge_x = 0.005)  +
scale_color_discrete(guide=FALSE) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + coord_cartesian(xlim=c(-0.15,0.16)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none")


#install.packages('ggrepel')
library(ggrepel)
ppr <- ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point(size=3) + geom_text_repel()  + scale_color_discrete(guide=FALSE)  +
coord_cartesian(xlim=c(-0.05,0.05),ylim=c(-0.05,0)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none")
ppf<-plot_grid(ppl, ppr, labels = c("A", "B"))
save_plot("~/tmp/bb_plot.pdf",ppf,ncol = 2,height=5)


if(FALSE){
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
}
