## here we compute a lookup that converts an observation of a single case genotype
## into an expected log(or)

library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)

or.threshold <- 2
nsims <- 1e6
n.sample <- 2500
or.prior <- 0.01

support.dir<-'/scratch/ob219/as_basis/support_tab'
## gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
## ref_af_file<-file.path(support.dir,'as_basis_snps.tab')
## ld_file<-file.path(support.dir,'all.1cM.tab')
## m_file<-file.path(support.dir,'as_basis_manifest.tab')
## ## dir where all preprocessed gwas files are.
## ## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
## basis.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir)
## basis.DT[,abs.or:=exp(abs(log(or)))]
## # what is the probability that a given snp has a or > threshold
## or.prior <- mean(sapply(split(basis.DT,basis.DT$trait),function(sDT){
##      mean(sDT$abs.or>or.threshold)
## }))
fname <- sprintf("lor_posterior_%g_%g_%g_%g.tab",or.threshold,or.prior,n.sample,nsims)
if(!file.exists(file.path(support.dir,fname))){
  af <- seq(0.01,0.99,by=0.01)
  library(parallel)
  options(mc.cores=6)
  results <- data.table(cbind(f=af,do.call("rbind",mclapply(af, function(f){
    message(sprintf("Processing %f",f))
    lor_f(f,n.sample,nsims,or.threshold,or.prior)
  }))))
  write.table(results,file=file.path(support.dir,fname),sep="\t",row.names=FALSE,quote=FALSE)
}else{
  results <- fread(file.path(support.dir,fname),header=TRUE)
}



library(ggplot2)
library(cowplot)
colnames(results) <- make.unique(colnames(results))
m <- reshape2::melt(results,"f")
ppf<-ggplot(m[grep("%", m$variable,invert=TRUE),],aes(x=f,y=value,col=variable,group=variable)) +
  geom_point() + labs(x="Allele Frequency",y="Posterior log(OR)") + background_grid(major = "xy", minor = "none") +
  scale_color_manual(name="Genotype",labels=c("0/0","0/1 or 1/0","1/1"),values=c("11"='firebrick',"01"='black',"00"='dodgerblue')) + theme(legend.position=c(0.1,0.8))
ppf
save_plot("~/tmp/poster_lor_plot.pdf",ppf)
