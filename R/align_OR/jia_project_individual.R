library(data.table)
library(annotSnpStats)
library(parallel)
library(cupcake)
support.dir<-'/scratch/ob219/as_basis/support_tab'
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
ref_af_file<-'/scratch/ob219/as_basis/support_tab/as_basis_snp_support.tab'
out_dir <- '/scratch/ob219/as_basis/jia_ind_analysis/'

or.threshold <- 10
nsims <- 1e6
n.sample <- 2500
or.prior <- 0.05
support.dir <- '/scratch/ob219/as_basis/support_tab'
lor.lu.file <- sprintf("lor_posterior_%g_%g_%g_%g.tab",or.threshold,or.prior,n.sample,nsims) %>%
    file.path(support.dir,.)


out_dir <- sprintf("ind_proj_%g_%g_%g_%g",or.threshold,or.prior,n.sample,nsims) %>%
    file.path(out_dir,.)
dir.create(out_dir, showWarnings = FALSE)


comp<-function(cv){
  tmp<-cv
  a<-c('A','G','C','T')
  b<-c('T','C','G','A')
  idx<-split(1:length(cv),cv)
  for(i in seq_along(a)){
    cv[idx[[a[i]]]]<-b[i]
  }
  cv
}
# flip snps
bsnps <- fread(ref_af_file)
setkey(bsnps,pid)
fls<-list.files(path='/scratch/wallace/JIA-2017-data/',pattern="*.RData",full.names=TRUE)
mclapply(fls,function(f){
  message(sprintf("Processing %s",f))
  load(f)
  sdt <- data.table(snps(G))[,c('lid','pid'):=list(1:.N,paste(chromosome,position,sep=':'))]
  setkey(sdt,pid)
  sdt <- sdt[bsnps][!(is.na(ref_a1.af) | is.na(chromosome)),]
  sdt <- sdt[,.(pid,allele.1,allele.2,ref_a1,ref_a2,ref_a1.af,lid)] %>%
    setnames(.,c('pid','a1','a2','ref_a1','ref_a2','ref_a1.af','lid'))
  keep <- sdt[sdt$pid %in% bsnps$pid,]$lid
  G <- G[,keep]
  ## check that thing are the same
  anno.snps<-data.table(snps(G))
  if(!identical(with(anno.snps,paste(chromosome,position,sep=':')),sdt$pid))
  stop(sprintf("Something went wrong for %s after subseting snp annotations don't line up",f))
  ## work out how we need to alter things so that we line up with reference snps
  ## note here we expect that outcome will be wrt to allele 2 (as snpStats)
  controls <- which(samples(G)$phenotype==0)
  sdt[,af.wrt.a2:=col.summary(G[controls,])[,"RAF"]]
  sdt[,flag:='unprocessed']
  sdt[a1==ref_a1 & a2==ref_a2,flag:='match']
  sdt[a1==ref_a2 & a2==ref_a1,flag:='flip']
  sdt[a1==comp(ref_a1) & a2==comp(ref_a2),flag:='match_rc']
  sdt[a1==comp(ref_a2) & a2==comp(ref_a1),flag:='flip_rc']
  ## Chris' code for assigning lor to individuals
  dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
  setnames(dat,c('f','lor.00','lor.01','lor.11'))
  lor <- as.matrix(dat[,-1])  # matrix, rows indexed by 100*AF, columns by gt 00/01/11

  ## next get genotypes for cases and project. From snpStats 1,2,3
  cases.idx<-which(samples(G)$phenotype==1)
  sm <- as(G,'SnpMatrix')
  sm<-matrix(sm[cases.idx,],nrow=length(cases.idx),ncol=ncol(sm))
  gt.lor <- lapply(seq_along(sdt$af.wrt.a2),function(j){
    if(sdt$flag[j] %in% c('flip','flip_rc')){
      g <- abs(3-as.numeric(sm[,j])) + 1
      i <- round((sdt$af.wrt.a2[j])*100)
    }else{
      g <- as.numeric(sm[,j])
      i <- round((1-sdt$af.wrt.a2[j])*100)
    }
    lor[i,g]
  }) %>% do.call("rbind",.)
  obj <- list(snps=sdt,proj.lor=gt.lor,samples=samples(G)[cases.idx,])
  fname <- gsub("annotsnpstats-([^.]+)\\.RData","chr\\1.RDS",basename(f))
  saveRDS(obj,file=file.path(out_dir,fname))
},mc.cores=8)


if(FALSE){
  library(data.table)
  library(ggplot2)
  dat<-readRDS(file.path(out_dir,"chr22.RDS"))
  sman<-dat$snps
  sman[,mean.proj.lor:=rowMeans(dat$proj.lor)]

  ## get the computed GWAS summary stats
  jia.DT<-fread('/home/ob219/scratch/as_basis/gwas_stats/input_files//jiacc_unpub.tab')[,.(pid,beta=log(or))]
  setkey(jia.DT,pid)
  setkey(sman,pid)
  sman <- jia.DT[sman]
  ## this line 'fixes things' by aligning the summary(OR) to the basis to fit the projections
  #sman[flip==TRUE,beta:=-beta] # to check
  ggplot(sman,aes(x=beta,y=mean.proj.lor,color=flag)) + geom_point() + geom_smooth(method="lm") +
  geom_abline(intercept=0,slope=1,color='red',lty=2) + theme_bw() + xlab("Case/Control Beta") +
  ylab("Mean(Beta.proj)")
  ggplot(sman,aes(x=beta,y=mean.proj.lor)) + geom_point() + geom_smooth(method="lm") +
  geom_abline(intercept=0,slope=1,color='red',lty=2) + theme_bw() + xlab("Case/Control Beta") +
  ylab("Mean(Beta.proj)")
}
