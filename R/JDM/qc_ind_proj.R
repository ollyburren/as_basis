## code to process the JDM summary statistics

## first load in summary data

library(data.table)
library(magrittr)

deak <- fread("/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/deakin/basicassoc_maf1_291016.csv")
deak[,pid:=paste(CHR,BP,sep=':')]
setkey(deak,pid)
## read in basis

support.dir<-'/home/ob219/rds/hpc-work/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support_feb_2018_w_ms.tab')
bsnps <- fread(ref_af_file)
setkey(bsnps,pid)

inc <- deak[bsnps][!is.na(STAT),]
## read in individual level data

all.ind.files <- list.files("/home/ob219/rds/hpc-work/as_basis/jdm_ind_analysis/",pattern="*.RDS",full.names=TRUE)

im.DT <- lapply(all.ind.files,function(f){
  dat <- readRDS(f)
  s <- dat$snps
  s[,mlor:=rowMeans(dat$proj.lor)]
  s
}) %>% rbindlist



ind <- im.DT[,.(pid,mlor)]
setkey(ind,pid)
inc[,sum.lor:=log(OR)]
## mark flips
inc[,flip:=A1!=ref_a2 & A1==ref_a1]
## fix flips !
summ <- inc[flip==TRUE,sum.lor:=sum.lor * -1][,.(pid,ref_a1.af,sum.lor,flip)]
setkey(summ,pid)
m<-ind[summ]

m[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
m[,mafl:=cut(m$maf,c(0,0.01,0.1,0.25,0.5))]

library(cowplot)

ggplot(m,aes(x=mlor,y=sum.lor)) + geom_point() + facet_wrap(~mafl) + geom_smooth(method="lm") + xlab("Individual log(POR)") + ylab("Summary Stats log(OR)")
