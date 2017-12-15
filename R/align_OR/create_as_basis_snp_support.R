library(data.table)
library(magrittr)
library(parallel)
library(GenomicRanges)

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


PROCESS_DIR<-'/scratch/ob219/as_basis/gwas_stats/processed/'
manifest <- fread('/scratch/ob219/as_basis/support_tab/as_manifest_december.tab')[basis_trait==1,]
ref <- fread('/scratch/ob219/as_basis/support_tab/EUR_1kg_support.tab')
MISSING_FILE_DEFAULT<-'/scratch/ob219/as_basis/support_tab/additional_missing_from_basis.txt'
ld_file<-'/scratch/ob219/as_basis/support_tab/all.1cM.tab'
maf.threshold <- 0.01

setkey(ref,pid)
setnames(ref,c('pid','ref_a1','ref_a2','ref_a1.af'))
## now attempt to align

all.snps <- rbindlist(mclapply(manifest$disease,function(d){
  message(d)
  ifile <- file.path(PROCESS_DIR,sprintf("%s.tab",d))
  ## remove things that are not SNPs when we load in
  input <- fread(ifile)[nchar(a1)==1 & nchar(a2)==1,]
  ## remove duplicates
  input <- input[!duplicated(pid),]
  setkey(input,pid)
  ## finally check that alleles line up
  out <- input[ref][!is.na(ref_a1.af),]
  out[,flag:='unprocessed']
  out[a1==ref_a1 & a2==ref_a2,flag:='match']
  out[a1==ref_a2 & a2==ref_a1,flag:='flip']
  out[a1==comp(ref_a1) & a2==comp(ref_a2),flag:='match_rc']
  out[a1==comp(ref_a2) & a2==comp(ref_a1),flag:='flip_rc']
  out[flag!='unprocessed',]
},mc.cores=8))

setkey(all.snps,pid)
summary <- all.snps[,list(s.count=.N),by=pid][s.count== nrow(manifest),]
## ok tally with reference file and save
out.ref <- ref[pid %in% summary$pid,]
## next remove MHC ! We define this broadly as chr6:20e6:40e6
coords <- tstrsplit(out.ref$pid,':')
snps.gr <- GRanges(seqnames=Rle(coords[[1]]),ranges=IRanges(start=as.numeric(coords[[2]]),width=1L))
ol <- findOverlaps(snps.gr,GRanges(seqnames=Rle('6'),ranges=IRanges(start=20e6,end=40e6))) %>% as.matrix
out.ref <- out.ref[-ol[,1],]
snps.gr <- snps.gr[-ol[,1],]
## add ld.blocks (this use to be done by cupcake but can do here instead)
if(!file.exists(ld_file))
  stop(sprintf("Cannot find file %s",ld_file))
ld<-fread(ld_file)
ld.gr<-with(ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
ol<-as.matrix(findOverlaps(snps.gr,ld.gr))
out.ref[ol[,1],ld.block := ol[,2]]
## finally remove any SNPs that are in missing file
if(file.exists(MISSING_FILE_DEFAULT)){
  missing<-scan(MISSING_FILE_DEFAULT,character())
  out.ref <- out.ref[!pid %in% missing,]
}
## finally remove snps that have a reference allele frequency less than 10%

maf<-with(out.ref,ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af))
out.ref <- out.ref[maf>maf.threshold,]
write.table(out.ref,file='/scratch/ob219/as_basis/support_tab/as_basis_snp_support.tab',quote=FALSE,row.names=FALSE,sep="\t")
