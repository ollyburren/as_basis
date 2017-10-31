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

gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)


study1<-'jia_cc'
study2<-'jia_ERA'
study1.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,study1)
study2.DT<-get_gwas_data(m_file,ref_af_file,ld_file,gwas_data_dir,study2)
## project onto basis to get the distance on PC2
jia.DT<-rbind(study1.DT,study2.DT)
setkey(jia.DT,pid)
jia.mat.emp<-create_ds_matrix(jia.DT,shrink.DT,'emp')
jia.pred.emp <- predict(pc.emp,newdata=jia.mat.emp)
observed.distance.pc2<-abs(jia.pred.emp[2,2]-jia.pred.emp[1,2])

## next compute the stderr for each - we use the se_maf_emp from the main basis
study1.DT <- study1.DT[shrink.DT[,.(pid,emp_maf_se)]]
study2.DT <- study2.DT[shrink.DT[,.(pid,emp_maf_se)]]

## next compute variance of beta for each SNP

emp_se <- function(DT){
  DT[,emp_se := 1/sqrt(2) * 1/sqrt(n) * emp_maf_se]
}

emp_se(study1.DT)
emp_se(study2.DT)
## here is some code to create snpMatrix objects for all snps in a basis
if(FALSE){
  bcft <- '~/bin/bcftools-1.4/bcftools'
  vcf.dir<-'/home/ob219/scratch/DATA/1kgenome/VCF/EUR/by.chr.phase1/'
  out.dir <- '/home/ob219/scratch/as_basis/snpStats/basis_1kg/'
  s.DT <- split(study1.DT[,.(chr,position)],study1.DT$chr)
  for(chr in names(s.DT)){
      message(sprintf("Processing %s",chr))
      tmp<-s.DT[[chr]][order(as.numeric(chr),as.numeric(position)),]
      tmp.file <- tempfile(pattern = "vcf_tmp", tmpdir = tempdir())
      write.table(tmp,file=tmp.file,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      #vcf.file <- file.path(vcf.dir,sprintf("ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz",chr))
      vcf.file <- file.path(vcf.dir,sprintf("chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz",chr))
      obj <- vcf2snpmatrix(vcf.file,bcft,tmp.file,TRUE)
      unlink(tmp.file)
      out.file <- file.path(out.dir,sprintf("%s_1kg.RData",chr))
      save(obj,file=out.file)
    }
}

#simulations
snpstats.dir <- '/home/ob219/scratch/as_basis/snpStats/basis_1kg/'
study1.DT[,or:=1]
study2.DT[,or:=1]
s1.sim <- simulate_study(study1.DT,snpstats.dir,shrink_beta=FALSE,n_sims=200,quiet=FALSE)
s2.sim <- simulate_study(study2.DT,snpstats.dir,shrink_beta=FALSE,n_sims=200,quiet=FALSE)
## think this should be done under the null where beta = 0 or or = 1

basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
## project on biobank to see if we can recreate Chris' figure.

sim.DT<-rbind(s1.sim,s2.sim)
setkey(sim.DT,pid)
sim.mat.emp<-create_ds_matrix(sim.DT,shrink.DT,'emp')
pred.emp <- predict(pc.emp,newdata=sim.mat.emp)
tnames<-row.names(pred.emp)
emp.DT<-data.table(pred.emp)
emp.DT[,trait:=tnames]
emp.DT <- melt(emp.DT,id.vars='trait')
emp.DT <- emp.DT[,c('bt','type','sim'):=tstrsplit(trait,'_')]
## compute delta distance for PC2
pc2.DT <- dcast(emp.DT[variable=='PC2',],sim~variable+type)[,delta:=PC2_ERA - PC2_cc]
emp.p<-sum(pc2.DT$delta>observed.distance.pc2)/nrow(pc2.DT)
## looks as if emp.p < 1/200 therefore this distance is significant.


## compute actual distance
emp<-rbind(pc.emp$x,pred.emp)
emp<-rbind(emp,jia.pred.emp)
ml<-list(
  CD = 'bb_CD',
  CEL = 'bb_CEL',
  MS = 'bb_MS',
  RA = 'bb_RA',
  SLE = 'bb_SLE',
  T1D = 'bb_T1D',
  UC = 'bb_UC',
  jia_cc = 'jia_ERA'
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
emp[grepl("ERA_V",trait),c('compare','trait'):=list('ERA','')]
emp[grepl("cc_V",trait),c('compare','trait'):=list('CC','')]
ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('JIA_ERA and cc comparison wrt to null')

# show distribution of distances

plot(density(pc2.DT$delta))
ggplot(pc2.DT,aes(x=delta)) + geom_histogram()
