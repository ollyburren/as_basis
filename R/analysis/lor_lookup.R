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

if(!file.exists(file.path(support.dir,'lor_posterior.tab'))){
  af <- seq(0.01,0.99,by=0.01)
  library(parallel)
  options(mc.cores=6)
  results <- data.table(cbind(f=af,do.call("rbind",mclapply(af, function(f){
    message(sprintf("Processing %f",f))
    lor_f(f,n.sample,nsims,or.threshold,or.prior)
  }))))
  write.table(results,file=file.path(support.dir,'lor_posterior.tab'),sep="\t",row.names=FALSE,quote=FALSE)
}else{
  results <- fread(file.path(support.dir,'lor_posterior.tab'))
}



library(ggplot2)
colnames(results) <- make.unique(colnames(results))
m <- reshape2::melt(results,"f")
ggplot(m[grep("%", m$variable,invert=TRUE),],
       aes(x=f,y=value,col=variable,group=variable)) + geom_point() + geom_smooth() + labs(x="AF",y="log OR")
