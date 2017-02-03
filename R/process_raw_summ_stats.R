## BEGAN: 18/01/2017
## PURPOSE: Process summary statistics for AS basis project.
library(data.table)
## simple format is id,chr,risk.allele,alt.allele,OR,p.value
raw.dir<-'/scratch/ob219/as_basis/gwas_stats/raw/'
process.dir<-'/scratch/ob219/as_basis/gwas_stats/processed/'
out.cols<-c('id','chr','position','risk.allele','other.allele','or','p.val')

flip<-function(DT){
	flip.idx<-which(DT$or<1)
	tmp<-DT[flip.idx,]$a1
	DT[flip.idx,]$a1<-DT[flip.idx,]$a2
	DT[flip.idx,]$a2<-tmp
	DT[flip.idx,]$or<-signif(1/DT[flip.idx,]$or,digits=5)
	setnames(DT,c('a1','a2'),c('risk.allele','other.allele'))
	DT
}

formatOut<-function(DT){
	DT<-setcolorder(DT,out.cols)
	DT[order(DT$chr,DT$position),]
}

processIMB<-function(DT){
	setnames(DT,c('id','chr','position','p.val','or','low','high','alleles'))
	tmp<-strsplit(DT$alleles,">")
	tmp<-matrix(unlist(tmp),ncol=2,byrow=TRUE)
	DT<-cbind(DT,tmp)
	setnames(DT,c('V1','V2'),c('a2','a1'))
	DT$or<-as.numeric(DT$or)
	DT$alleles<-NULL
	DT<-flip(DT)
	DT<-DT[,c('id','chr','position','p.val','or','other.allele','risk.allele'),with=FALSE]
	formatOut(DT)
}

processIBD<-function(DT){
	DT$or<-exp(DT$Effect)
	DT<-DT[,c('MarkerName','Allele1','Allele2','P.value','or'),with=FALSE]
	DT$chr<-sub("([^:]+):.*","\\1",DT$MarkerName)
	DT$position<-sub("[^:]+:([0-9]+).*","\\1",DT$MarkerName)
	DT$position<-as.numeric(sub("[^:]+:([0-9]+).*","\\1",DT$MarkerName))
	##convert alleles to upper case
	for(n in c('Allele1','Allele2'))
	        DT[[n]]<-toupper(DT[[n]])
	setnames(DT,c('id','a1','a2','p.val','or','chr','position'))
	DT<-flip(DT)
	DT<-formatOut(DT)
}

## T2D
t2d<-fread(paste0(raw.dir,'DIAGRAMv3.2012DEC17.txt'))
# build appears to be 36 therefore move to 37
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
alleles.loc <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37,t2d$SNP,ifnotfound="drop")
t2d$b37.position <- as.numeric(NA)
setkey(t2d,'SNP')
t2d[alleles.loc$RefSNP_id,]$b37.position <- start(alleles.loc)
# assume OR is wrt to risk allele
if(any(t2d$OR<1)){
	message("Help OR is not always wrt to risk allele")
}
## could check alleles
t2d.out<-subset(t2d,!is.na(b37.position))
t2d.out<-t2d.out[,c('SNP','CHROMOSOME','RISK_ALLELE','OTHER_ALLELE','P_VALUE','OR','b37.position'),with=FALSE]
setnames(t2d.out,c('id','chr','risk.allele','other.allele','p.val','or','position'))
t2d.out<-formatOut(t2d.out)
save(t2d.out,file=file.path(process.dir,'T2D.RData'))
rm(list=c('t2d.out','t2d','alleles.loc'))

##RA
ra<-fread(paste0(raw.dir,'RA_GWASmeta_European_v2.txt'))
setnames(ra,c('id','chr','position','a1','a2','or','hi','low','p.val'))
ra.out<-ra[,c('id','chr','position','a1','a2','or','p.val'),with=FALSE]
ra.out<-flip(ra.out)
ra.out<-formatOut(ra.out)
save(ra.out,file=file.path(process.dir,'RA.RData'))

##SCZ
scz<-fread(paste0(raw.dir,'scz2.snp.results.txt'))
setnames(scz,c('chr','id','a1','a2','position','info','or','se','p.val','ngt'))
## assume that the or is wrt to a1
scz.out<-scz[,c('id','chr','position','a1','a2','or','p.val'),with=FALSE]
scz.out$chr<-as.numeric(substr(scz.out$chr,4,nchar(scz.out$chr)))
scz.out<-flip(scz.out)
save(scz.out,file=file.path(process.dir,'SCZ.RData'))
rm(list=c('scz','scz.out','ra','ra.out'))

##CEL
cel<-fread(paste0(raw.dir,'hg19_gwas_cel_dubois_4_19_1.tab'))
## all immunobase files can probably be processed the same way
cel.out<-processIMB(cel)
save(cel.out,file=file.path(process.dir,'CEL.RData'))

##SLE
sle<-fread(paste0(raw.dir,'hg19_gwas_sle_bentham_4_20_0.tab'))
sle.out<-processIMB(sle)
sle.out<-subset(sle.out,!is.na(sle.out$or))
save(sle.out,file=file.path(process.dir,'SLE.RData'))

##MS
ms<-fread(paste0(raw.dir,'hg19_gwas_ms_imsgc_4_19_2.tab'))
ms.out<-processIMB(ms)
ms.out<-subset(ms.out,!is.na(ms.out$or))
save(ms.out,file=file.path(process.dir,'MS.RData'))

##UC
uc<-fread(paste0(raw.dir,'uc_build37_45975_20161107.txt'))
uc.out<-processIBD(uc)
save(uc.out,file=file.path(process.dir,'UC.RData'))

##CD
cd<-fread(paste0(raw.dir,'cd_build37_40266_20161107.txt'))
cd.out<-processIBD(cd)
save(cd.out,file=file.path(process.dir,'CD.RData'))

##IBD
ibd<-fread(paste0(raw.dir,'ibd_build37_59957_20161107.txt'))
ibd.out<-processIBD(ibd)
save(ibd.out,file=file.path(process.dir,'IBD.RData'))

##PSCG
psc<-fread(paste0(raw.dir,'ipscsg2016.result.combined.full.with_header.txt'))
psc.out<-psc[,c('#chr','SNP','pos','allele_0','allele_1','or','p'),with=FALSE]
setnames(psc.out,c('chr','id','position','risk.allele','other.allele','or','p.val'))
setcolorder(psc.out,out.cols)
save(psc.out,file=file.path(process.dir,'PSC.RData'))

##T1D
t1d.dir<-'/scratch/wallace/gwas-summary-stats/T1DGC+WTCCC-1KG-imputed/'
t1d.files<-list.files(path=t1d.dir,patter='snptest-meta.*',full.names=TRUE)


t1d<-rbindlist(lapply(seq_along(t1d.files),function(i){
	tf<-t1d.files[i]
	chr<-sub(".*chr([^\\.]+)\\.gz","\\1",basename(tf))
	message(chr)
	out<-fread(paste('zcat',tf))
	out$chr<-chr
	out
}))
t1d$chromosome<-t1d$chr
t1d$chr<-NULL

## function to process T1D and get different cohorts based on colnames 

processT1D<-function(t1d,beta.cname,se.cname){
	cnames<-c('rsid','chromosome','position','alleleA','alleleB',beta.cname,se.cname)
	tmp<-t1d[,cnames,with=FALSE]
	setnames(tmp,c(cnames[1:5],'beta','se'))
	tmp$or<-exp(tmp$beta)
	tmp$Z<-abs(tmp$beta/tmp$se)
	tmp$p<-pnorm(-abs(tmp$Z))
	tmp.out<-tmp[,c('rsid','chromosome','position','alleleA','alleleB','or','p'),with=FALSE]
	setnames(tmp.out,c('id','chr','position','a1','a2','or','p.val'))
	tmp.out<-flip(tmp.out)
	formatOut(tmp.out)
}

t1d.meta<-processT1D(t1d,'beta.meta','se.meta')
save(t1d.meta,file=file.path(process.dir,'T1D_meta.RData'))
t1d.1<-processT1D(t1d,'beta.1','se.1')
save(t1d.1,file=file.path(process.dir,'T1D_1.RData'))
t1d.2<-processT1D(t1d,'beta.2','se.2')
save(t1d.2,file=file.path(process.dir,'T1D_2.RData'))


## meta
#t1d$or<-exp(t1d$beta.meta)
#t1d$Z<-abs(t1d$beta.meta/t1d$se.meta)
#t1d$p<-pnorm(-abs(t1d$Z))
#t1d.out<-t1d[,c('rsid','chromosome','position','alleleA','alleleB','or','p'),with=FALSE]
#setnames(t1d.out,c('id','chr','position','a1','a2','or','p.val'))
#setcolorder(t1d.out,out.cols)
#save(t1d.out,file=file.path(process.dir,'T1D.RData'))
