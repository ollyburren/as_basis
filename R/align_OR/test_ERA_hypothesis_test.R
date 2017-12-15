library(data.table)

simulation.dir <- '/home/ob219/scratch/as_basis/jia_summ_analysis'
basis.cache.file <- '~/tmp/as_basis_cache.RDS'
trait <- 'ERA'
fp <- sprintf("jia_%s*",trait)
sim.files <- list.files(path=simulation.dir,pattern=fp,full.names=TRUE)
all.sims <- do.call('rbind',(lapply(sim.files,readRDS)))
## read in the actual basis
basis <- readRDS(basis.cache.file)
t.test(all.sims[,"PC2"],mu=basis$basis$x["control","PC2"])

## for ERA jia ?

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
jia_traits<-fread(m_file)[grep('jia',trait),]$trait
jia.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,jia_traits)
jia.mat.emp<-create_ds_matrix(jia.DT,shrink.DT,'emp')
pred.emp <- predict(pc.emp,newdata=jia.mat.emp)

t.test(all.sims[,"PC2"],mu=pred.emp["jia_ERA","PC2"])

sum(all.sims[,"PC2"]>pred.emp["jia_ERA","PC2"])/nrow(all.sims)
