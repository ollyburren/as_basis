library(data.table)
## this is a script to take summary stats from Astle et al and format so we can project onto AI basis

DATA.DIR <- '/scratch/wallace/gwas-summary-stats/blood/ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2016-12-12/ukbiobank'
support.dir<-'/scratch/ob219/as_basis/support_tab'
ref_af_file<-file.path(support.dir,'as_basis_snp_support.tab')
OUT.DIR <- '/home/ob219/scratch/as_basis/gwas_stats/raw_blood_new/'
thou.file <-'/home/ob219/scratch/as_basis/1KG_support/all_EUR_0.01_support.RData'
if(FALSE){
  ## first we load in and filter so that we have just the SNPs we need for the basis
  ## quickest way is to use tabix

  ## first use ref_af_file to create a temp file
  (load(thou.file))
  #DT.s <- fread(ref_af_file)[,c('chr','start'):=tstrsplit(pid,':')]
  DT.s <- all.eur[,.(chr,position)]
  DT.s <- DT.s[order(as.numeric(chr),as.numeric(position))]
  target.file <- file.path('/scratch/ob219/as_basis/support_tab','tabix_support.tab')
  write.table(DT.s[,.(chr,position)],file=target.file,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)


  TABIX_BIN<-'~/bin/htslib/tabix -R %s %s > %s '
  bfiles <- list.files(path=DATA.DIR,pattern="*.gz$",full.names=TRUE)
  cmds<-lapply(bfiles,function(f){
    of <- file.path(OUT.DIR,gsub("\\.gz","",basename(f)))
    sprintf(TABIX_BIN,target.file,f,of)
  })
  write(unlist(cmds),file="~/tmp/run_blood_tabix.txt")
  #run on the queue using qlines.R
}
tab.head <- c('VARIANT','ID','CHR','BP','REF',
'ALT','ALT_MINOR','DIRECTION','EFFECT',
'SE','P','MLOG10P','ALT_FREQ','MA_FREQ')


library(devtools)
library(mvtnorm)
devtools::install_github("ollyburren/cupcake")
library(cupcake)
library(magrittr)

IN.DIR <- '/home/ob219/scratch/as_basis/gwas_stats/raw_blood_new/'
bfiles <- list.files(path=IN.DIR,pattern="*.tsv$",full.names=TRUE)
OUT.DIR <- '/home/ob219/scratch/as_basis/gwas_stats/processed_new/'

for(f in bfiles){
  message(sprinf("Processing %s",basename(f)))
  DT <- fread(f)
  setnames(DT,tab.head)
  ## extract trait and sample size
  fname <- basename(f)
  trait <- gsub("(.*)\\_build37\\_[0-9]+\\_20161212.tsv","\\1",fname)
  ss <- as.numeric(gsub("(.*)\\_build37\\_([0-9]+)\\_20161212.tsv","\\2",fname))
  smaf <- DT$MA_FREQ
  sbeta <- DT$EFFECT
  sbeta.se <- DT$SE
  res.DT <- convertBetaToOR(N=ss,b=sbeta,seb=sbeta.se,m=smaf) %>% do.call('cbind',.) %>% data.table

  ## with processed datasets we need to convert into the correct format.
  #pid,
  #a1
  #a2
  #or
  #p.val

  setnames(DT,'P','P.cont')
  all.DT <- cbind(DT,res.DT)
  all.DT[,pid:=paste(CHR,BP,sep=':')]
  out <- all.DT[,.(pid,REF,ALT,signif(OR,digits=4),signif(P,digits=4))]
  setnames(out,c('pid','a1','a2','or','p.val'))
  of <- file.path(OUT.DIR,sprintf("astle_%s.tab",trait))
  write.table(out,file=of,quote=FALSE,sep="\t",row.names=FALSE)
  message(sprinf("Wrote %s",of))
}
