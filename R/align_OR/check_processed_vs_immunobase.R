# Get curations from ImmunoBase for diseases and assign effect alleles.

IMMUNOBASE_URL_TEMPLATE <- 'https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-%s-assoc_variantsTAB'
DISEASE2MESH <- list(
  SLE = 'Lupus Erythematosus, Systemic',
  CRO = 'Crohn Disease',
  CEL = 'Celiac Disease',
  MS = 'Multiple Sclerosis',
  PSC = 'Cholangitis, Primary Sclerosing',
  PBC = 'Liver Cirrhosis, Biliary',
  UC = 'Colitis, Ulcerative',
  T1D = 'Diabetes Mellitus, Type 1',
  RA = 'Arthritis, Rheumatoid'
)
PROCESS_DIR <- '/scratch/ob219/as_basis/gwas_stats/processed/'
col.names<-c('pid','a1','a2','or','p.value')
library(data.table)
library(magrittr)

get_curation_DT<-function(disease){
  url <- sprintf(IMMUNOBASE_URL_TEMPLATE,disease)
  DT <- fread(sprintf("curl -s %s",url)) %>% setnames(.,make.names(names(.)))
  DT <- DT[Disease==DISEASE2MESH[[disease]],]
  DT[,c('a1','a2'):=tstrsplit(Alleles,'>')]
  DT[,chr:=gsub("([^:]+):.*","\\1",DT$Coords.)]
  DT[chr=='X',chr:='23']
  DT[,pid:=paste(chr,Position,sep=':')]
  DT <- DT[,.(pid,a1,a2,Odds.Ratio,P.Value)] %>% setnames(.,col.names)
  DT[,direction:=sign(log(or))]
  setkey(DT,pid)
  DT
}

compare_or <- function(target,reference,flip=FALSE){
  setkey(target,pid)
  target[,direction:=sign(log(or))]
  tmp <- target[reference]
  message(nrow(tmp))
  if(!flip){
    match <- subset(tmp,direction==i.direction & i.a1==a1 & i.a2==a2)
    message(sprintf("Matched %d",nrow(match)))
    subset(tmp,direction!=i.direction & i.a1==a1 & i.a2==a2)
  }else{
    match <- subset(tmp,direction==i.direction & i.a2==a1 & i.a1==a2)
    message(sprintf("Matched %d",nrow(match)))
    subset(tmp,direction!=i.direction & i.a2==a1 & i.a1==a2)
  }
}

#sle

sle <- fread(file.path(PROCESS_DIR,'sle_bennett.tab'))
sle.cur <- get_curation_DT('SLE')
compare_or(sle,sle.cur)
## looks as if errors in ImmunoBase OR flipped for the following.

#cd
cd <- fread(file.path(PROCESS_DIR,'cd_delaange.tab'))
cd.cur <- get_curation_DT('CRO')
compare_or(cd,cd.cur)

#cel
cel <- fread(file.path(PROCESS_DIR,'cel_dubois.tab'))
cel.cur <- get_curation_DT('CEL')
compare_or(cel,cel.cur)

#ms
ms <- fread(file.path(PROCESS_DIR,'ms_imsgc.tab'))
ms.cur <- get_curation_DT('MS')
compare_or(ms,ms.cur,flip=TRUE)

#psc
psc <- fread(file.path(PROCESS_DIR,'psc_ji.tab'))
psc.cur <- get_curation_DT('PSC')
compare_or(psc,psc.cur)

#pbc
pbc <- fread(file.path(PROCESS_DIR,'pbc_cordell.tab'))
pbc.cur <- get_curation_DT('PBC')
compare_or(pbc,pbc.cur)

#uc
uc <- fread(file.path(PROCESS_DIR,'uc_delaange.tab'))
uc.cur <- get_curation_DT('UC')
compare_or(uc,uc.cur)

#t1d
t1d <- fread(file.path(PROCESS_DIR,'t1d_cooper.tab'))
t1d.cur <- get_curation_DT('T1D')
compare_or(t1d,t1d.cur)

##ra (for completeness)
ra <- fread(file.path(PROCESS_DIR,'ra_okada.tab'))
ra.cur <- get_curation_DT('RA')
compare_or(ra,ra.cur,flip=TRUE)

## next we compute posterior probabilities and flag regions where credible SNP interval is 1

## load in manifest main is on gdocs https://docs.google.com/spreadsheets/d/1v2R6ehdqanM3WJjKhXLNQB_2M2e7EvvD712aAj_2U1E/edit?usp=sharing
library(cupcake)
ref_1kg <- fread('/home/ob219/scratch/as_basis/support_tab/EUR_1kg_support.tab')
manifest <- fread('/home/ob219/scratch/as_basis/support_tab/as_manifest_december.tab')
manifest[,N:=cases+controls]
ld <- fread('/home/ob219/scratch/as_basis/support_tab/all.1cM.tab')
ld.gr <- GRanges(seqnames=Rle(ld$chr),ranges=IRanges(start=ld$start,end=ld$end))
setkey(ref_1kg,pid)
setnames(ref_1kg,c('pid','ref.a1','ref.a2','a1.af'))

## locate suspicious regions

sus_reg <- function(d,thresh.ppi=0.9,thresh.pval=1e-5){
  disease.DT <- manifest[label==d,]
  ifile <- file.path(PROCESS_DIR,sprintf("%s.tab",disease.DT$disease))
  message(sprintf("Reading %s",ifile))
  tmp <- fread(ifile)
  setkey(cd,pid)
  tmp<-ref_1kg[tmp][!is.na(ref.a2),]
  coords<-tstrsplit(tmp$pid,':')
  ol <- GRanges(seqnames=Rle(coords[[1]]),ranges=IRanges(start=as.numeric(coords[[2]]),width=1L)) %>% findOverlaps(.,ld.gr) %>% as.matrix
  tmp[ol[,1],ld.block:=ol[,2]]
  tmp[,maf:=ifelse(a1.af>0.5,1-a1.af,a1.af)]
  ##next compute wakefield

  tmp[,ppi:=wakefield_pp(p.value,maf,disease.DT$N,disease.DT$cases/disease.DT$N,1e-4,0.2),by=ld.block]
  tmp <- tmp[,.SD[order(ppi,decreasing=TRUE)],by=ld.block]
  tmp[,c.sum:=cumsum(ppi),by=ld.block]
  ## define some threshold for ppi say 0.9 and see howmany snps are in the credible set

  ci.sets <- tmp[,list(ci.set=length(which(c.sum<thresh.ppi & .N>1))),by=ld.block]
  ## get a list of those that look suspicious
  tmp[ld.block %in% ci.sets[ci.set==1,]$ld.block,list(sig.p=length(which(p.value<thresh.pval)),no=.N),by=ld.block]
}

foo<-lapply(manifest[basis_trait==1,]$label,sus_reg)

## can't see anything majorly off.
