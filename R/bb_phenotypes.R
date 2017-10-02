library(data.table)

scratch.dir<-'/home/ob219/scratch/as_basis'
bb_pheno_url <- 'https://www.dropbox.com/s/oe5q85454vhc3hi/phenosummary_final_11898_18597.tsv?dl=0'
bb_pheno_file <- file.path(scratch.dir,'bb','phenosummary_final_11898_18597.tsv')
if(!file.exists(bb_pheno_file)){
  system(sprintf('wget %s -O %s',bb_pheno_url,bb_pheno_file))
}
pheno <- fread(bb_pheno_file)

## use self reported instead. The ICD classifications are incomplete meaning that the controls are polluted
bb_codes<-c('20002_1226','20002_1225','20002_1111','20002_1452','20002_1222','20002_1261','20002_1313','20002_1462','20002_1459','20002_1463','20002_1464','20002_1381','20002_1456')

 ## get these
curated<-subset(pheno,Field.code %in% bb_codes)[,.(Field.code,Field,N.cases,N.controls)]
curated$Field<-sub('^[^:]+:[ ](.*)$','\\1',curated$Field)
curated$cat<-TRUE


all.self.report.cc<-subset(pheno,!is.na(N.cases) & grepl("Non-cancer illness code, self-reported",Field) & !Field.code %in% bb_codes & Field != 'unclassifiable')[,.(Field.code,Field,N.cases,N.controls)]
all.self.report.cc$Field<-sub('^[^:]+:[ ](.*)$','\\1',all.self.report.cc$Field)
all.self.report.cc$cat<-FALSE

## plot to see the distro of background case controls with ai case controls
library(ggplot2)
ggplot(rbind(all.self.report.cc[all.self.report.cc$N.cases>min(curated$N.cases),],curated),aes(x=cat,y=N.cases)) + geom_boxplot()

## we decided to take the top 10 (excluding unclassifiable)

background<-head(all.self.report.cc[order(all.self.report.cc$N.cases,decreasing=TRUE),],n=10)

## get the list of links from here (this of course could change)
#https://docs.google.com/spreadsheets/d/1b3oGI2lUt57BcuHttWaZotQcI0-mBRPyZihz87Ms_No/edit#gid=1209628142

bb.links<-fread(file.path(scratch.dir,'bb','bb_gwas_linklist.20170915.csv'))
setnames(bb.links,make.names(names(bb.links)))

## get the data for the traits that we want

bb.get<-subset(bb.links,Phenotype.code %in% c(background$Field.code,curated$Field.code))

## get and convert to RDS format
# db.url = drop box address
# fc = field code used to compose out file
# out.dir = directory where RDS will be written
convert2RDS<-function(db.url,fc,desc,out.dir=scratch.dir){
  ## get the data as a tmp file
  tfile<-tempfile(fileext='.gz')
  system(sprintf('wget %s -O %s',db.url,tfile))
  tmp<-fread(sprintf('gunzip -c %s',tfile))
  ncols<-length(names(tmp))
  nrows<-nrow(tmp)
  fname<-paste(fc,make.names(desc),sep=':')
  ofile<-file.path(out.dir,'bb','summary_stats',paste(fname,'RDS',sep='.'))
  saveRDS(tmp,file=ofile)
  file.remove(tfile)
  return(data.table(ofile=ofile,n.cols=ncols,n.rows=nrows))
}


outcome<-lapply(1:nrow(bb.get),function(i){
  df<-bb.get[i,]
  message(sprintf("Processing %s",df$Description))
  with(df,convert2RDS(db.url=Dropbox.address,fc=Phenotype.code,desc=Description))
})















## this code is old but may be useful at a later stage
get_icd<-function(){
    ## these were prepared from table 1 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892249/
    icd<-readRDS("/home/ob219/scratch/as_basis/support/ai_icd.RDS")[,.(disease,icd10)]

    ## in BB don't seem to have classifications after the .

    icd$icd10simp<-sub("([^.]*)\\..*","\\1",icd$icd10)
    icd<-icd[!duplicated(icd$icd10simp),]

    ## get ai traits by ICD

    all.res<-lapply(1:nrow(icd),function(i){
      grepterm<-paste('main ICD10:',icd[i,]$icd10simp)
      print(grepterm)
      pheno[grep(grepterm,pheno$Field,ignore.case=TRUE),]
    })
    names(all.res)<-icd$disease
    phe.det<-rbindlist(all.res)

    ## get self reported

    tl<-tolower(sub("Diagnoses - main ICD10: [A-Z][0-9][0-9]*","",phe.det$Field))

    ## some manual sorting here as codes are a bit of a mess

    keyword<-c('hypothyroidism','hyperthyroidism','type 1 diabetes','multiple sclerosis','crohns disease',
    'guillain-barre',)

    subset(pheno,!is.na(N.cases))
}
