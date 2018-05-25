    library(data.table)
    library(annotSnpStats)
    library(parallel)
    library(cupcake)
    support.dir<-'/scratch/ob219/as_basis/support_tab'
    gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
    ref_af_file<-file.path(support.dir,'as_basis_snps_with_af_fixed.tab')
    out_dir <- '/scratch/ob219/as_basis/jia_ind_analysis/ind_proj_aligned_fixed'
    # flip snps
    flips<-readRDS("/home/ob219/scratch/as_basis/tmp/support_file_flips.RDS")
    bsnps <- fread(ref_af_file)
    fls<-list.files(path='/scratch/wallace/JIA-2017-data/',pattern="*.RData",full.names=TRUE)
    mclapply(fls,function(f){
      message(sprintf("Processing %s",f))
      load(f)
      sdt <- data.table(snps(G))
      keep <- which(paste(sdt$chromosome,sdt$position,sep=':') %in% bsnps$pid)
      G <- G[,keep]
      #test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(G),snp.data=as(G,"SnpMatrix"))
      test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(G),snp.data=as(G,"SnpMatrix"))
      DF<-do.call('rbind',test)
      DF<- apply(DF,2,unlist)
      sman<-data.table(snps(G))
      sman[,c('beta','Var.beta'):=list(as.numeric(DF[,"beta"]),as.numeric(DF[,"Var.beta"]))]
      controls <- which(samples(G)$phenotype==0)
      sman[,af.wrt.a2:=col.summary(G[controls,])[,"RAF"]]
      ## work out which are flipped wrt to the basis
      sman[,pid:=paste(chromosome,position,sep=':')]
      align.DT <- sman[,.(chromosome,position,allele.1,allele.2,pid)]
      setnames(align.DT,c('chr','position','a1','a2','pid'))
      if(FALSE){
        sman[,flip:=FALSE]
        sman[flip_allele(align.DT,bsnps),flip:=TRUE]
      }else{
        sman[,flip:=pid %in% flips]
      }
      ## Chris' code for assigning lor to individuals
      lor.lu.file <- file.path(support.dir,'lor_posterior.tab')

      #bsnps <- fread(ref_af_file)
      #bsnps <- bsnps[sample(1:nrow(bsnps),10000),]
      dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
      setnames(dat,c('f','lor.00','lor.01','lor.11'))
      lor <- as.matrix(dat[,-1])  # matrix, rows indexed by 100*AF, columns by gt 00/01/11

      ## next get genotypes for cases and project. From snpStats 1,2,3
      cases.idx<-which(samples(G)$phenotype==1)
      library(magrittr)
      sm <- as(G,'SnpMatrix')
      sm<-matrix(sm[cases.idx,],nrow=length(cases.idx),ncol=ncol(sm))
      gt.lor <- lapply(seq_along(sman$af.wrt.a2),function(j){
        if(sman$flip[j]){
          g <- abs(3-as.numeric(sm[,j])) + 1
          i <- round((sman$af.wrt.a2[j])*100)
        }else{
          g <- as.numeric(sm[,j])
          i <- round((1-sman$af.wrt.a2[j])*100)
        }
        lor[i,g]
      }) %>% do.call("rbind",.)
      obj <- list(snps=sman,proj.lor=gt.lor,samples=samples(G)[cases.idx,])
      fname <- gsub("annotsnpstats-([^.]+)\\.RData","chr\\1.RDS",basename(f))
      saveRDS(obj,file=file.path(out_dir,fname))
    },mc.cores=8)


if(FALSE){
  library(data.table)
  library(ggplot2)
  out_dir <- '/scratch/ob219/as_basis/jia_ind_analysis/ind_proj_aligned_fixed'
  out_dir <- '/home/ob219/rds/hpc-work/as_basis/jia_ind_analysis/ind_proj_aligned'
  dat<-readRDS(file.path(out_dir,"chr22.RDS"))
  sman<-dat$snps
  sman[,mean.proj.lor:=rowMeans(dat$proj.lor)]
  ## this line 'fixes things' by aligning the summary(OR) to the basis to fit the projections
  #sman[flip==TRUE,beta:=-beta] # to check
  ggplot(sman,aes(x=beta,y=mean.proj.lor,color=flip)) + geom_point() + geom_smooth(method="lm") +
  geom_abline(intercept=0,slope=1,color='red',lty=2) + theme_bw() + xlab("Case/Control Beta") +
  ylab("Mean(Beta.proj)")
}
