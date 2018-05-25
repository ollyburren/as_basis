library(data.table)
library(magrittr)

## first off create a manifest of all snps that have been impute but using cut on the gen files
info.files <- list.files(path="/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed",pattern="imputed-wtccc-.*.gen_filt_info$",full.names=TRUE)

## read in the current basis support file.

support.DT<-fread("/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab")
support.DT[,c('chr','position'):=tstrsplit(pid,':')]

by.chr <- split(support.DT$pid,support.DT$chr)

miss<-list()

OUT.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed/tmp'

gtool_bin <- "/home/ob219/bin/gtool/gtool"

# gtool -S --g example/example.gen --s example/example.sample --og example/out.gen --os example/out.sample --inclusion example/rs_id.txt

cmd_pat <- "%s -S --g %s  --og %s --inclusion %s"

DAT.DIR <- "/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/t1d-barrett-imputed"
gstuff <- lapply(info.files,function(f){
  chr<-gsub("imputed-wtccc-([^\\.]+)\\.gen_filt_info","\\1",basename(f))
  DT<-fread(sprintf("cut -d' ' -f2,3 %s",f))[,c('pid','line'):=list(paste(chr,position,sep=':'),1:.N)]
  fDT <- subset(DT,pid %in% by.chr[[chr]])
  if(nrow(fDT) != length(by.chr[[chr]])){
    tmp<-by.chr[[chr]]
    message(sprintf("Missing %d",length(by.chr[[chr]])-nrow(fDT)))
    #miss[[chr]]<-tmp[!tmp %in% fDT$pid]
  }
  fDT
  fout <- file.path(OUT.DIR,gsub("gen_filt_info","txt",basename(f)))
  write(fDT$rs_id,file=fout)
  ## next formulate gtools command
  in.file <- file.path(DAT.DIR,sprintf("imputed-wtccc-%s.gen.gz",chr))
  out.file <- file.path(OUT.DIR,sprintf("t1d_bootstrap_imputed-wtccc-%s.gen.gz",chr))
  fDT[,cmd:=sprintf(cmd_pat,gtool_bin,in.file,out.file,fout)]
  fDT
}) %>% rbindlist

## get a list of missing snps

new.support.DT <- support.DT[pid %in% gstuff$pid,]
write.table(new.support.DT,file="/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_t1d_bootstrap.tab",sep="\t",quote=FALSE,row.names=FALSE)

comm<-unique(gstuff$cmd)
write(comm,file=file.path(OUT.DIR,"run.txt"))
