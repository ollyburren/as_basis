library(data.table)
library(GenomicRanges)



REMOVE_MHC <- TRUE
INTERSECT_FILE <- '~/tmp/intersect_w_ms.txt'
# name,chr,position,a1,a2,a2,f,pid
REFERENCE_FILE <- '/home/ob219/scratch/as_basis/1KG_support/all_EUR_0.01_support.RData'
LD_BLOCK_FILE <- '/home/ob219/scratch/as_basis/support_tab/all.1cM.tab'
SUPPORT_FILE <- '/scratch/ob219/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
all.eur <- get(load(REFERENCE_FILE))
keep.pid <- scan(INTERSECT_FILE,character())
ref <- all.eur[pid %in% keep.pid, ]
ref[,.(pid,ref_a1=a1,ref_a2=a2,read_a1.af=1-a2.f)]

## remove MHC
if(REMOVE_MHC){
  snps.gr <- with(ref,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,width=1L),pid=pid))
  mhc.gr <- GRanges(seqnames=Rle('6'),ranges=IRanges(start=20e6,end=40e6))
  mhc.pid <- subsetByOverlaps(snps.gr,mhc.gr)$pid
  ref <- ref[!pid %in% mhc.pid,]
}

## add in LD assignments
snps.gr <- with(ref,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,width=1L),pid=pid))
if(!file.exists(LD_BLOCK_FILE))
  stop(sprintf("Cannot find file %s",LD_BLOCK_FILE))
ld<-fread(LD_BLOCK_FILE)
ld.gr<-with(ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
ol<-as.matrix(findOverlaps(snps.gr,ld.gr))
ref[ol[,1],ld.block := ol[,2]]
ref <- ref[order(chr,position),.(pid,ref_a1=a1,ref_a2=a2,ref_a1.af=1-a2.f,ld.block)]
write.table(ref,file=SUPPORT_FILE,quote=FALSE,row.names=FALSE,sep="\t")
