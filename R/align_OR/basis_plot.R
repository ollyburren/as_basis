library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)
## import GWAS data for basis
## support files

SHRINKAGE_METHOD<-'ws_emp'

support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support.tab')
#ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-file.path(support.dir,'as_manifest_december.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'

support.dir<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
# reference allele frequencies
ref_af_file<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
#ld_file<-file.path(support.dir,'all.1cM.tab')
#m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv'
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/rds/user/ob219/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'

basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
emp<-pc.emp$x
vexp <- summary(pc.emp)[['importance']][2,]
PC1.var<-signif(vexp["PC1"]*100,digits=3)
PC2.var<-signif(vexp["PC2"]*100,digits=3)
M <- cbind(as.data.table(emp),trait=rownames(emp))
scp <- cbind(data.table(vexp),pcs=factor(names(vexp),levels=names(vexp)))
scp[,cs:=cumsum(vexp)]
scp$group=1

## nice plot
library("cowplot")
library("ggrepel")

## do a scree plot

ppl <- ggplot(scp,aes(x=pcs,y=vexp,group=group)) + geom_point() + geom_line() + ylab("Variance Explained") + xlab("Principal Components")

ppr<-ggplot(M,aes(x=PC1,y=PC2,label=trait)) + geom_point(size=3) + geom_text_repel() + # hjust = 0, nudge_x = 0.005)  +
scale_color_discrete(guide=FALSE) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + coord_cartesian(xlim=c(-0.15,0.16)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none")


ppf<-plot_grid(ppl, ppr, labels = c("A", "B"))
save_plot("~/tmp/basis_plot.pdf",ppf,ncol = 2,base_height=5)
