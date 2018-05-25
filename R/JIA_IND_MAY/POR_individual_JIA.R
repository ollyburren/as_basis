 library(data.table)
 library(annotSnpStats)
 library(parallel)
 library(cupcake)


 ## code to create the Posterior Odds for JIA dataset

 gwas_data_dir <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
 ref_af_file<-file.path("/home/ob219/rds/rds-cew54-wallace-share/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab")
 #basis.snps <- fread(ref_af_file)[,pid:=paste(chr,position,sep=':')]
 basis.snps <- fread(ref_af_file)[,c('chr','position'):=tstrsplit(pid,":")]
 setkey(basis.snps,pid)

 #support.dir<-'/scratch/ob219/as_basis/support_tab'
 #gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
 #ref_af_file<-file.path(support.dir,'as_basis_snps_with_af_fixed.tab')
 out_dir <- '/home/ob219/rds/hpc-work/as_basis/jia_ind_analysis'
 # flip snps
 bsnps <- fread(ref_af_file)
 fls<-list.files(path='/home/ob219/rds/hpc-work/as_basis/JIA_basis_annotSnpStats/',pattern="*.RData",full.names=TRUE)

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

 getAFLU <- function(af){
   if(between(af,1,99))
     return(af)
   if(af==0)
     return(1)
   if(af==100)
     return(99)
 }


lapply(fls,function(f){
   message(sprintf("Processing %s",f))
   load(f)
   sdt <- data.table(snps(G))
   sdt[,pid:=paste(chromosome,position,sep=':')]
   ## remove non SNPs
   rm.idx <- which(nchar(as.character(sdt$a1))!=1)
   if(length(rm.idx)>0){
     sdt[rm.idx,pid:=paste(pid,"remove",sep='_')]
     colnames(G) <- sdt$pid
     G <- G[,-rm.idx]
     sdt <- sdt[-rm.idx,]
   }
   # duplicates
   rm.idx <- which(duplicated(sdt$pid))
   if(length(rm.idx)>0){
     sdt[rm.idx,pid:=paste(pid,"dup",sep='_')]
     colnames(G) <- sdt$pid
     G <- G[,-rm.idx]
     sdt <- sdt[-rm.idx,]
   }
   ## need to set alleles to be correct
   alleles(G) <- c('a0','a1')
   colnames(G) <- sdt$pid
   keep <- which(paste(sdt$chromosome,sdt$position,sep=':') %in% bsnps$pid)
   G <- G[,keep]
   #test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(G),snp.data=as(G,"SnpMatrix"))
   #test=snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(G),snp.data=as(G,"SnpMatrix"))
   #DF<-do.call('rbind',test)
   #DF<- apply(DF,2,unlist)
   #sman<-data.table(snps(G))
   #sman[,c('beta','Var.beta'):=list(as.numeric(DF[,"beta"]),as.numeric(DF[,"Var.beta"]))]
   #controls <- which(samples(G)$phenotype==0)
   #sman[,af.wrt.a2:=col.summary(G[controls,])[,"RAF"]]
   ## work out which are flipped wrt to the basis
   sman <- data.table(snps(G))[,ord:=1:.N]
   sman[,rs_id:=as.character(rs_id)]
   sman[,pid:=paste(chromosome,position,sep=':')]
   setkey(bsnps,pid)
   setkey(sman,pid)
   sman<-bsnps[sman]
   setnames(sman,c('a0','a1'),c('a1','a2'))
   sman[,flag:='unprocessed']
   sman[a1==ref_a1 & a2==ref_a2,flag:='match']
   sman[a1==ref_a2 & a2==ref_a1,flag:='flip']
   sman[a1==comp(ref_a1) & a2==comp(ref_a2),flag:='match_rc']
   sman[a1==comp(ref_a2) & a2==comp(ref_a1),flag:='flip_rc']
   ## here ref_af2 is wrong !!
   sman[,tmp:=af.wrt.a2]
   sman[,af.wrt.a2:=1-ref_a1.af]



   sman <- sman[order(ord),]
   #align.DT <- sman[order(ord),.(chromosome,position,a1,a2,pid,flag)]
   #setnames(align.DT,c('chr','position','a1','a2','pid','flag'))

   ## Chris' code for assigning lor to individuals
   support.dir<-'/home/ob219/rds/hpc-work/as_basis/support_tab'
   lor.lu.file <- file.path(support.dir,'lor_posterior.tab')

   #bsnps <- fread(ref_af_file)
   #bsnps <- bsnps[sample(1:nrow(bsnps),10000),]
   dat <- fread(lor.lu.file,skip=1L)[,.(V1,V2,V3,V4)]
   setnames(dat,c('f','lor.00','lor.01','lor.11'))
   lor <- as.matrix(dat[,-1])  # matrix, rows indexed by 100*AF, columns by gt 00/01/11

   ## next get genotypes for cases and project. From snpStats 1,2,3
   #cases.idx<-which(samples(G)$phenotype==1)
   library(magrittr)
   sm <- as(G,'SnpMatrix')
   sm<-matrix(sm,nrow=nrow(sm),ncol=ncol(sm))


   # this deals with imputed genotypes
   gt.lor <- lapply(seq_along(sman$af.wrt.a2),function(j){
     if(sman$flag[j] %in% c('flip','flip_rc')){
       #i <- getAFLU(round((1-sman$af.wrt.a2[j])*100))
       i <- getAFLU(round((1-sman$af.wrt.a2[j])*100))
       as.vector(pp(sm[,j]) %*% rev(lor[i,]))
     }else{
       #i <- getAFLU(round((1-sman$af.wrt.a2[j])*100))
       i<-getAFLU(round((1-sman$af.wrt.a2[j])*100))
       as.vector(pp(sm[,j]) %*% lor[i,])
     }
   }) %>% do.call("rbind",.)
   # gt.lor <- lapply(seq_along(sman$af.wrt.a2),function(j){
   #   message(j)
   #   if(sman$flip[j]){
   #     g <- abs(3-as.numeric(sm[,j])) + 1
   #     i <- round((sman$af.wrt.a2[j])*100)
   #   }else{
   #     g <- as.numeric(sm[,j])
   #     i <- round((1-sman$af.wrt.a2[j])*100)
   #   }
   #   lor[i,g]
   # }) %>% do.call("rbind",.)
   obj <- list(snps=sman,proj.lor=gt.lor,samples=samples(G))
   fname <- gsub("annotsnpstats-([^.]+)\\.RData","chr\\1.RDS",basename(f))
   saveRDS(obj,file=file.path(out_dir,fname))
 })
