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
	## note that we switch a1 and a2 as alleles are wrt to a2
	setnames(DT,c('id','a2','a1','p.val','or','chr','position'))
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

##PBC
pbc<-fread(paste0(raw.dir,'hg19_gwas_pbc_cordell_4_20_0.tab'))
## a bunch (7219) of pbc don't appear to have proper alleles from some reason
## perhaps go back and fix by taking raw data and remapping to b37 for time being
## filter out

pbc<-pbc[pbc[[8]] != '?/?',]
pbc.out<-processIMB(pbc)
pbc.out<-subset(pbc.out,!is.na(pbc.out$or))
save(pbc.out,file=file.path(process.dir,'PBC.RData'))

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
#setnames(psc.out,c('chr','id','position','risk.allele','other.allele','or','p.val'))
setnames(psc.out,c('chr','id','position','other.allele','risk.allele','or','p.val'))
setcolorder(psc.out,out.cols)
save(psc.out,file=file.path(process.dir,'PSC.RData'))

##T1D - Compared to Nick's this is wrong all OR are halved - DONT USE
#t1d.dir<-'/scratch/wallace/gwas-summary-stats/T1DGC+WTCCC-1KG-imputed/'
#t1d.files<-list.files(path=t1d.dir,patter='snptest-meta.*',full.names=TRUE)


#t1d<-rbindlist(lapply(seq_along(t1d.files),function(i){
#	tf<-t1d.files[i]
#	chr<-sub(".*chr([^\\.]+)\\.gz","\\1",basename(tf))
#	message(chr)
#	out<-fread(paste('zcat',tf))
#	out$chr<-chr
#	out
#}))
#t1d$chromosome<-t1d$chr
#t1d$chr<-NULL

## function to process T1D and get different cohorts based on colnames

#processT1D<-function(t1d,beta.cname,se.cname){
#	cnames<-c('rsid','chromosome','position','alleleA','alleleB',beta.cname,se.cname)
#	tmp<-t1d[,cnames,with=FALSE]
#	setnames(tmp,c(cnames[1:5],'beta','se'))
#	tmp$or<-exp(tmp$beta)
#	tmp$Z<-abs(tmp$beta/tmp$se)
#	tmp$p<-pnorm(-abs(tmp$Z))
#	tmp.out<-tmp[,c('rsid','chromosome','position','alleleA','alleleB','or','p'),with=FALSE]
#	setnames(tmp.out,c('id','chr','position','a1','a2','or','p.val'))
#	tmp.out<-flip(tmp.out)
#	formatOut(tmp.out)
#}
#
#t1d.meta<-processT1D(t1d,'beta.meta','se.meta')
#save(t1d.meta,file=file.path(process.dir,'T1D_meta.RData'))
#t1d.1<-processT1D(t1d,'beta.1','se.1')
#save(t1d.1,file=file.path(process.dir,'T1D_1.RData'))
#t1d.2<-processT1D(t1d,'beta.2','se.2')
#save(t1d.2,file=file.path(process.dir,'T1D_2.RData'))

##T1D these are derived from Nick Coopers 1KG Imputed GWAS
t1d.file<-'/scratch/ob219/ncooper_bioarxiv_2017/raw/meta_all.txt'
t1d.out<-fread(t1d.file)
## filter as only want variants that are in all basis
setkeyv(t1d.out,c('chromosome','position'))
setnames(t1d.out,'id','tid')
t1d.out$chromosome<-as.numeric(t1d.out$chromosome)
## we use asthma to filter only to variants found in other studies
if(!exists("asthma.out"))
	load(file.path(process.dir,'asthma.RData'))
asthma.out$chr<-as.numeric(asthma.out$chr)
#tmp<-t1d.out[asthma.out] - this doesn't work not sure why
tmp<-merge(t1d.out,asthma.out,by.x=c('chromosome','position'),by.y=c('chr','position'))[,names(t1d.out),with=FALSE]

processT1D<-function(t1d,or.cname,se.cname,p.cname){
	cnames<-c('rsid','chromosome','position','a0','a1',or.cname,se.cname,p.cname)
	tmp<-t1d[,cnames,with=FALSE]
	setnames(tmp,c(cnames[1:5],'or','se','p.val'))
	#tmp$p.val<-2*pnorm(abs(log(tmp$or)/tmp$se),lower.tail=FALSE)
	tmp.out<-tmp[!is.na(tmp$p.val),c('rsid','chromosome','position','a0','a1','or','p.val'),with=FALSE]
	## note that OR is wrt allele2 so flip
	setnames(tmp.out,c('id','chr','position','a2','a1','or','p.val'))
	tmp.out<-flip(tmp.out)
	formatOut(tmp.out)
}

aff.t1d<-processT1D(tmp,'aff.OR','aff.se','aff.pvalue')
save(aff.t1d,file=file.path(process.dir,'aff.t1d.RData'))
ill.t1d<-processT1D(tmp,'ill.OR','ill.se','ill.pvalue')
save(ill.t1d,file=file.path(process.dir,'ill.t1d.RData'))
meta.t1d<-processT1D(tmp,'OR.meta','se.meta','p.meta')
save(meta.t1d,file=file.path(process.dir,'meta.t1d.RData'))

## JIA
jia.dir<-'/scratch/wallace/gwas-summary-stats/jia-2017-unpublished/'
jia.files<-list.files(path=jia.dir,pattern='*.gz',full.names=TRUE)

jia<-rbindlist(lapply(jia.files,function(f){
	message(paste("Processing",f))
	tmp<-fread(paste('zcat',f),skip=1L)
	setnames(tmp,c('id','chr','snp.name','cM','position','a1','a2','other.id','ref','alt','af','beta.sys','se.sys','beta.nosys','se.nosys'))
	tmp[,c('id','chr','position','a1','a2','beta.sys','se.sys','beta.nosys','se.nosys'),with=FALSE]
}))


processJIA<-function(jia,beta.cname,se.cname){
	cnames<-c('id','chr','position','a1','a2',beta.cname,se.cname)
        tmp<-jia[,cnames,with=FALSE]
        setnames(tmp,c(cnames[1:5],'beta','se'))
        tmp$or<-exp(tmp$beta)
        tmp$Z<-abs(tmp$beta/tmp$se)
        tmp$p<-pnorm(-abs(tmp$Z))
        tmp.out<-tmp[,c('id','chr','position','a1','a2','or','p'),with=FALSE]
	## note that or is wrt allele2 sp flip
        setnames(tmp.out,c('id','chr','position','a2','a1','or','p.val'))
        tmp.out<-flip(tmp.out)
        formatOut(tmp.out)
}

jia.sys<-processJIA(jia,'beta.sys','se.sys')
save(jia.sys,file=file.path(process.dir,'JIA_sys.RData'))
jia.nosys<-processJIA(jia,'beta.nosys','se.nosys')
save(jia.nosys,file=file.path(process.dir,'JIA_nosys.RData'))

## Asthma
asthma.dir<-'/scratch/wallace/gwas-summary-stats/asthma/'
a.file<-file.path(asthma.dir,'gabriel.csv.gz')
as.DT<-fread(sprintf("zcat %s",a.file))
as.DT<-as.DT[,c('Chr','rs','position','Allele_1','Allele_2','OR_ran','ORl_ran','ORu_ran','P_ran'),with=FALSE]
# cannot do anything where we don't have the odds ratio
as.DT<-as.DT[!is.na(as.DT$OR_ran),]
## double check things
#as.DT[,`:=`(uSE=abs(log(ORu_ran)-log(OR_ran))/1.96,lSE=abs(log(ORl_ran)-log(OR_ran))/1.96),]
#as.DT$uCalcP<-with(as.DT,2*pnorm(abs(log(OR_ran)/uSE),lower.tail=FALSE))
#as.DT$lCalcP<-with(as.DT,2*pnorm(abs(log(OR_ran)/lSE),lower.tail=FALSE))
#cor(as.DT$uCalcP,as.DT$P_ran)
#cor(as.DT$lCalcP,as.DT$P_ran)
## remap to build 37 -- what is best way liftover or using an annotation library ?
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
snps<-SNPlocs.Hsapiens.dbSNP144.GRCh37
lu<-snpsById(snps,as.DT$rs,ifnotfound='drop')
lu<-data.table(as.data.frame(lu))
setkey(lu,'RefSNP_id')
setkey(as.DT,'rs')
as.DT<-as.DT[lu]
as.DT<-as.DT[,c('rs','seqnames','pos','Allele_1','Allele_2','OR_ran','P_ran'),with=FALSE]
setnames(as.DT,c('id','chr','position','a1','a2','or','p.val'))
as.DT$chr<-sub("^ch","",as.DT$chr)
as.DT<-flip(as.DT)
asthma.out<-formatOut(as.DT)
save(asthma.out,file=file.path(process.dir,'asthma.RData'))

## blood traits - these come from Astle et al. but they have 30M variants we only require a few of these.
blood.data.dir<-'/scratch/wallace/gwas-summary-stats/blood/ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2016-12-12/ukbiobank'
bfiles<-c('eo_build37_172275_20161212.tsv.gz','lymph_build37_171643_20161212.tsv.gz','myeloid_wbc_build37_169219_20161212.tsv.gz','wbc_build37_172435_20161212.tsv.gz')
## our basis can only use SNPs in common therefore filter datasets (we can use asthma above)
asthma.out$chr<-as.integer(asthma.out$chr)
setkeyv(asthma.out,c('chr','position'))
all.blood<-lapply(bfiles,function(f){
	message(sprintf("Processing %s",f))
	tmp<-fread(sprintf("zcat %s",file.path(blood.data.dir,f)))
	setkeyv(tmp,c('CHR','BP'))
	tmp<-tmp[asthma.out]
	tmp<-tmp[,names(tmp),with=FALSE]
	tmp<-tmp[!is.na(tmp$P),c('ID','CHR','BP','REF','ALT','EFFECT','P'),with=FALSE]
	## EFFECT is wrt ALT therefore need to swap. Also Effect is the log(OR) so need to compute exp(EFFECT)
	setnames(tmp,c('id','chr','position','a2','a1','or','p.val'))
	tmp$or<-exp(tmp$or)
	tmp<-flip(tmp)
	tmp<-formatOut(tmp)
})
names(all.blood)<-c('eosinophil','lymphocyte','myeloid','wbc')
for(i in seq_along(all.blood)){
	t<-names(all.blood)[i]
	varname<-paste(t,'out',sep='.')
	ofile<-file.path(process.dir,paste(t,'RData',sep='.'))
	assign(varname,all.blood[[i]])
	save(list=varname,file=ofile)
}

## I computed the GWAS summary statistics for JIA many subtypes

jia.data.dir<-'/home/ob219/scratch/jia/by.trait/'
jia.fs<-list.files(path=jia.data.dir,pattern='*.RDS',full.names=TRUE)
all.jia<-lapply(jia.fs,function(f){
	tmp<-readRDS(f)[,.(rsid,chr,position,a1,a2,beta,p.val)]
	## for SNP matrix all OR are wrt to a2
	tmp$or<-exp(tmp$beta)
	tmp.out<-tmp[,c('rsid','chr','position','a1','a2','or','p.val'),with=FALSE]
## note that or is wrt allele2 sp flip
	setnames(tmp.out,c('id','chr','position','a2','a1','or','p.val'))
	tmp.out<-flip(tmp.out)
	formatOut(tmp.out)
})
names(all.jia)<-paste('jia',basename(gsub("\\.RDS","",jia.fs)),sep='_')
for(i in seq_along(all.jia)){
	t<-names(all.jia)[i]
	varname<-paste(t,'out',sep='.')
	ofile<-file.path(process.dir,paste(t,'RData',sep='.'))
	assign(varname,all.jia[[i]])
	save(list=varname,file=ofile)
}

## Neale Lab Biobank summary statistics.
## We process these by a specific script process_bb_summary_stats.R

## Myositis

my.DT<-fread('/home/ob219/scratch/as_basis/gwas_stats/to_project/Myositis.txt')

## assume that these are all JDM
## I think that the data is from build36 so need to remap to 37
lu<-snpsById(snps,my.DT$SNP,ifnotfound='drop')
lu<-data.table(as.data.frame(lu))
setkey(lu,'RefSNP_id')
setkey(my.DT,'SNP')
my.DT<-my.DT[lu]
## attempt to reconstruct allele 2 from IUPAC codes notice there is some ambigutity.
cIUPAC<-list(
	R=c('A','G'),
	Y=c('C','T'),
	S=c('G','C'),
	W=c('A','T'),
	K=c('G','T'),
	M=c('A','C'),
	B=c('C','G','T'),
	D=c('A','G','T'),
	H=c('A','C','T'),
	N='N',
	V=c('A','C','G'))
)
