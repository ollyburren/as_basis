library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(optparse)

TEST<-TRUE
option_list = list(
        make_option(c("-s", "--snp_support_file"), type="character", default=NULL,
              help="SNP support file", metavar="character"),
        make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="Output file for basis", metavar="character"),
        make_option(c("-m", "--manifest_file"), type="character", default=NULL,
              help="Manifest file to use", metavar="character"),
        make_option(c("-g", "--gwas_data_dir"), type="character", default=NULL,
              help="Dir containing filtered, aligned GWAS summary stats",metavar="character")
        )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}else{
  args <- list(
      snp_support_file='/scratch/ob219/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab',
      outfile="~/tmp/testbasis.RDS",
      manifest_file='/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv',
      gwas_data_dir='/scratch/ob219/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
  )
}



basis.DT<-get_gwas_data(args$manifest_file,args$snp_support_file,args$gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

save(pc.emp,file=)

DT<-data.table(trait=rownames(pc.emp$x),PC1=pc.emp$x[,"PC1"],PC2=pc.emp$x[,"PC2"],basis.trait=T)
library(ggplot2)
ggplot(DT,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text()

## attempt to project on Lyons Data

lyons.traits <- fread(args$manifest_file)[grep("lyons|astle_eo$|astle_neut$",disease),]$trait
lyons.DT<-get_gwas_data(args$manifest_file,args$snp_support_file,args$gwas_data_dir,trait_list=lyons.traits)
lyons.mat.emp<-create_ds_matrix(lyons.DT,shrink.DT,'emp')
pred.emp <- predict(pc.emp,newdata=lyons.mat.emp)
DT.lyons <- data.table(trait=toupper(rownames(pred.emp)),PC1=pred.emp[,"PC1"],PC2=pred.emp[,"PC2"],basis.trait=F)
library(cowplot)
library(ggrepel)
ggplot(rbind(DT,DT.lyons),aes(x=PC1,y=PC2,label=trait,color=basis.trait)) + geom_point() +
geom_text_repel() + scale_color_manual(guide=FALSE,values=c('firebrick','dodgerblue'))
