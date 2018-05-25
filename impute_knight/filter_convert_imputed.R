## knight imputation data
library(data.table)

DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/expr/knight-imputed-genotypes/'
bsnps <- fread('/home/ob219/scratch/as_basis/support_tab/as_basis_snp_support.tab')
OUT.DIR <- file.path(DATA.DIR,'as_basis/IMPUTE')

## we want to filter the imputed knight data to get just the SNPs in the basis or
## as many of them as possible

pos <- do.call('rbind',strsplit(bsnps$pid,':'))
lu <- split(as.numeric(pos[,2]),pos[,1])

info.files <- list.files(path=DATA.DIR,pattern="*gen_filt_info",full.names=FALSE)
gen.files <- list.files(path=DATA.DIR,pattern="*gen_filt.gz",full.names=FALSE)
files <- data.table(info=info.files,gen=gen.files)

for(i in seq_along(1:nrow(files))){
  fin <- file.path(DATA.DIR,files$info[i])
  fgen <- file.path(DATA.DIR,files$gen[i])
  print(paste("Procesing",fin))
  tin <- fread(fin)
  chr <- gsub("imputed-knight-(.*).gen_filt_info","\\1",basename(fin))
  idx <- which(tin$position %in% lu[[chr]])
  tin <- tin[idx,]
  gen <- fread(paste('zcat < ',fgen))[idx,]
  gz1 <- gzfile(file.path(OUT.DIR,basename(fgen)), "w")
  options(scipen=999)
  write.table(gen,gz1,row.names=FALSE,sep=" ",col.names=FALSE,quote=FALSE)
  close(gz1)
  write.table(tin,file.path(OUT.DIR,basename(fin)),row.names=FALSE,sep=" ",col.names=TRUE,quote=FALSE)
}

## next read in and store as snpStats

library(annotSnpStats)
library(data.table)

convert2snpstats <- function(f.gt,info.filter=0.8){
    chr <- gsub("imputed-knight-([^\\.]+)\\..*","\\1",basename(f.gt))
    ss <- read.impute(f.gt)
    f.info <- gsub("\\.gz","_info",f.gt)
    ss.info <- fread(f.info)
    ss.info[,chr:=chr]
    ss.info[,c('rsid','position','allele.1','allele.2'):=tstrsplit(rs_id,':')]
    ss.info<-as.data.frame(ss.info[,.(rsid,allele.1,allele.2,chr,position,info)])
    ## add a filter on info to get rid of dodgy ones
    idx<-which(ss.info$info<info.filter)
    if(length(idx)>0){
      ss.info<-ss.info[-idx,]
      ss<-ss[,-idx]
    }
    rownames(ss.info)<-ss.info$rsid
    colnames(ss)<-ss.info$rsid
    samples <- data.frame(sample.id=rownames(ss))
    rownames(samples) <- rownames(ss)
    new("aSnpMatrix",
     .Data=ss,
     snps=ss.info,
     samples=samples )
}

AOUT.DIR <- file.path(DATA.DIR,'as_basis/annotSnpStats')
IN.DIR <- file.path(DATA.DIR,'as_basis/IMPUTE')
gen.files <- list.files(path=IN.DIR,pattern="*gen_filt.gz",full.names=TRUE)
lapply(gen.files,function(g){
  message(g)
  obj <- convert2snpstats(g)
  chr <- unique(snps(obj)$chr)
  fname<-file.path(AOUT.DIR,sprintf("chr%s.RDS",chr))
  saveRDS(obj,file=fname)
})

## next combine in one file
library(annotSnpStats)
ddir <- '/home/ob219/rds/rds-cew54-wallace-share/Data/expr/knight-imputed-genotypes/as_basis/annotSnpStats'
fs <- list.files(path=ddir,pattern="*.RDS",full.names=TRUE)
all.snps<-lapply(fs,readRDS)

all.gt <- do.call('cbind',lapply(all.snps,function(x) as(x,"SnpMatrix")))
all.info <- do.call('rbind',lapply(all.snps,snps))
names(all.info) <- c('rsid','A1','A2','chr','position','info')
combined <- new("aSnpMatrix",
 .Data=all.gt,
 snps=all.info,
 samples=samples(all.snps[[1]]))
saveRDS(combined,file="/home/ob219/rds/rds-cew54-wallace-share/Data/expr/knight-imputed-genotypes/as_basis/annotSnpStats/as_basis_knight_genotypes.RDS")

## first create snpSta
