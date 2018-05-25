library(data.table)
## need to get non astle and bb files and filter so that in 1k genome > 1%

bb_var_info<-readRDS('/home/ob219/scratch/as_basis/bb/variants_prefiltered_1kg_EUR_0.01.RDS')
bb_var_info[,c('chr','pos','a1','a2'):=tstrsplit(variant,':')]
bb_var_info[,pid:=paste(chr,pos,sep=':')]
OLD.DIR <- '/home/ob219/scratch/as_basis/gwas_stats/processed'
NEW.DIR <- '/home/ob219/scratch/as_basis/gwas_stats/processed_new'

o.f <- list.files(path=OLD.DIR,pattern="*.tab",full.names=TRUE)

## remove astle and biobank

to.process <- o.f[-grep("astle.tab|bb\\:",o.f)]

for(f in to.process){
  message(sprintf("Processing %s",f))
  DT <- fread(f)
  out <- DT[pid %in% bb_var_info$pid & !duplicated(pid),]
  write.table(out,file=file.path(NEW.DIR,basename(f)),quote=FALSE,sep="\t",row.names=FALSE)
}


## aav stuff

aav.dir <- '/rds/project/cew54/rds-cew54-wallace-share/Data/GWAS/egpa_aav/summary'

aav.files <- list.files(path=aav.dir,pattern='*.gwas',full.names=TRUE)
lapply(aav.files,function(f){
  DT <- fread(f)
  DT <- DT[order(CHR,BP),.(pid=paste(CHR,BP,sep=':'),a1=a0,a2=a1,or=signif(OR,digits=4),p.value=signif(P,digits=4))]
  out <- DT[pid %in% bb_var_info$pid & !duplicated(pid),]
  trait <- gsub("autosomes\\.(.*)\\.high-info\\.gwas","\\1",basename(f))
  fname <- sprintf("%s_lyons.tab",trait)
  write.table(out,file=file.path(NEW.DIR,fname),quote=FALSE,sep="\t",row.names=FALSE)
})
