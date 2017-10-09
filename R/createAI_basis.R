setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")
source("./wakefield.R")
library(magrittr)
sim.size<-100
pi_1<-1e-4
tmpfile<-"/home/ob219/scratch/as_basis/merged_data/AI_test_raw.RDS"

## these are the traits from which to build the basis
basis_traits<-c('CD','CEL','MS','PBC','PSC','RA','SLE','UC','meta.t1d')
#jia<-c('jia_cc','jia_EO','jia_ERA','JIA_nosys','jia_PO','jia_PsA','jia_RFneg','jia_RFpos','jia_sys')
# for ease of computation for simulation we will probably take forward gh_ss
# for computation without including quantitative traits

if(!file.exists(tmpfile)){
    DT<-getGWASData()
    ## remove quat traits
    #QT<-c('eosinophil','lymphocyte','myeloid','wbc')
    #RMT<-c(QT,'aff.t1d','ill.t1d')
    DT<-DT[DT$disease %in% basis_traits,]
    ## asthetically nicer if we move meta.t1d <- t1d
    DT[DT$disease=='meta.t1d',disease:='T1D']

    DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
    ## compute gamma - this can be computed in different ways
    ## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
    ##
    DT[,gh_ss:=gamma_hat_ss(Z,total)]
    DT[,gh_ss_pp:=gh_ss * pp]
    # examine what happens if we compute gh using maf ?
    DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
    DT[,gh_maf_pp:=gh_maf * pp]
    #DT<-DT[!DT$disease %in% jia,]
    DT[,lp0:=log(1-approx.bf.z2(Z,maf,cases+controls,cases/(cases+controls),pi_1)),by=c('disease','ld.block')]
    DT<-rbind(DT,createControl(DT))

    DT[DT$disease!='controls',q_i:=1-exp(sum(lp0)),by=id]
    emp<-mean(DT[,list(mean_qi=mean(q_i)),by=id]$mean_qi)
    ## prior odds
    po<-emp/(1-emp)
    ## note that emp us for h1 that beta != 0 therefore we need to take reciprocal as equation assumes pi_0 - see notes
    po<-1/po
    DT[,uABF:=po*(q_i/(1-q_i))]
    BIG<-quantile(DT[is.finite(DT$uABF),]$uABF,prob=0.9999)
    ## BIG is still rather large but at least tractable
    DT[is.infinite(uABF) | uABF > BIG, uABF:=BIG]
    DT[,pwi:=computePWI(uABF,emp),by=c('disease','ld.block')]
    DT[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]
    saveRDS(DT,file=tmpfile)
}else{
  DT<-readRDS(tmpfile)
}

## code to create matrix for any metric


createBasis<-function(DT,m){
  tmp<-melt(DT,id.vars=c('id','disease'),measure.vars = m)
  basis.matrix<-data.table::dcast(tmp,id~disease)
  basis.matrix<-basis.matrix[,2:ncol(basis.matrix),with=FALSE] %>% as.matrix(.) %>% t(.)
  prcomp(basis.matrix,center=TRUE,scale=FALSE)
}

## next we project on biobank data and see if it makes any sense whatsoever

all.DT<-getGWASData()
traits<-unique(all.DT$disease)
## get biobank data

bb.traits<-traits[grep('^bb',traits)]

bb.DT<-all.DT[all.DT$disease %in% bb.traits,]
bb.DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
## compute gamma - this can be computed in different ways
## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
##
bb.DT[,gh_ss:=gamma_hat_ss(Z,total)]
bb.DT[,gh_ss_pp:=gh_ss * pp]
# examine what happens if we compute gh using maf ?
bb.DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
bb.DT[,gh_maf_pp:=gh_maf * pp]
setkey(bb.DT,'id')
## need to add weightings
pwi<-unique(DT[,.(id,pwi)])
setkey(pwi,id)
bb.DT<-bb.DT[pwi]
bb.DT[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]


## metric/weighting to use
metric<-'gh_maf_pw'


## ok first we create basis
asb<-createBasis(DT,metric)

## for each biobank trait we project this onto the matrix

bb.preds<-lapply(seq_along(bb.traits),function(i){
  bbt<-bb.traits[i]
  message(sprintf("processing %s",bbt))
  proj<-subset(bb.DT,disease==bbt)[[metric]] %>% as.matrix(.) %>% t(.)
  pred<-data.table(predict(asb,proj))
  pred$trait <- bbt
  pred
})

bb.pc<-rbindlist(bb.preds)
bb.pc$basis<-FALSE
basis.pc <- data.table(asb$x)
basis.pc[,c('trait','basis'):=list(rownames(asb$x),TRUE)]
all.pc<-rbind(bb.pc,basis.pc)


## draw a scree plot

## draw one of these for each bb trait.

ml<-list(
  CD = 'bb:20002_1462:self_reported_crohns_disease',
  CEL = 'bb:20002_1456:self_reported_malabsorption_coeliac_disease',
  MS = 'bb:20002_1261:self_reported_multiple_sclerosis',
  RA = 'bb:20002_1464:self_reported_rheumatoid_arthritis',
  SLE = 'bb:20002_1381:self_reported_systemic_lupus_erythematosis_sle',
  T1D = 'bb:20002_1222:self_reported_type_1_diabetes',
  UC = 'bb:20002_1463:self_reported_ulcerative_colitis'
)

abpos<-c('PSC','RA','PBC','CEL','T1D','MS','SLE')
all.pc$biobank<-FALSE
all.pc[grep('^bb',all.pc$trait),biobank:=TRUE]
all.pc$ab.pos<-FALSE
all.pc[all.pc$trait %in% abpos | all.pc$trait %in% unlist(ml[abpos]),ab.pos:=TRUE]

m.pc<-melt(all.pc,id.vars=c('trait','ab.pos','biobank','basis'))

plotter<-function(i){
  basis_trait<-names(ml)[i]
  bb_trait<-ml[[i]]
  title<-paste(basis_trait,bb_trait)
  dat<-m.pc[which(m.pc$basis | m.pc$trait==bb_trait),]
  cf.idx<-which(dat$trait %in% c(basis_trait,bb_trait,'control'))
  dat$compare<-FALSE
  dat<-dat[cf.idx,compare:=TRUE]
  dat[dat$variable=='PC1',trait_label:=gsub('bb:[^:]+:self_reported_','',trait)]
  #ggplot(dat,aes(x=variable,y=value,lty=ab.pos,color=trait_label,group=trait,alpha=trait!='control')) + geom_point() + geom_line() + theme_bw()
  #ggplot(dat,aes(x=variable,y=value,label=trait_label,color=compare,group=trait,lty=ab.pos)) +
  ggplot(dat,aes(x=variable,y=value,label=trait_label,color=compare,group=trait,alpha=compare)) +
  #geom_point() + geom_line() + geom_text(angle=ifelse(dat$compare,90,0)) +
  geom_point() + geom_line() + geom_text() +
  scale_color_discrete(guide=FALSE)  + theme_bw() + ggtitle(title) + xlab('PC') + ylab('Loading')
}

pheno<-fread('/home/ob219/scratch/as_basis/gwas_stats/sample_counts_bb.csv')

sub.p<-subset(pheno,disease %in% c(names(ml),unlist(ml)))
sub.p<-sub.p[order(sub.p$cases,decreasing=TRUE),][,.(disease,cases)]
sub.p$disease<-gsub('bb:[^:]+:self_reported_','',sub.p$disease)

library(xtable)
xtable(pheno)

library(ggplot2)
pdf(file="~/tmp/bb_compare_as_basis.pdf")
lapply(seq_along(ml),function(i){
  plotter(i)
})
dev.off()

## define a plotter that allows us to hilight bb dataset on the basis

bb_plotter<-function(bb_trait){
  dat<-m.pc[which(m.pc$basis | m.pc$trait==bb_trait),]
  dat$hi<-FALSE
  dat<-dat[dat$trait %in% c(bb_trait,'control'),hi:=TRUE]
  dat[dat$variable=='PC1',trait_label:=gsub('bb:[^:]+:self_reported_','',trait)]
  ggplot(dat,aes(x=variable,y=value,label=trait_label,color=hi,group=trait,alpha=hi)) +
  geom_point() + geom_line() + geom_text() +
  scale_color_discrete(guide=FALSE)  + theme_bw() + ggtitle(bb_trait) + xlab('PC') + ylab('Loading')
}

pdf(file="~/tmp/bb_compare_all_as_basis.pdf")
lapply(unique(m.pc[!m.pc$basis,]$trait),bb_plotter)
dev.off()

## just check what happens if we reproject RA on to the basis.

ra.proj<-subset(DT,disease=='RA')[[metric]] %>% as.matrix(.) %>% t(.)
data.table(predict(asb,ra.proj))
all.pc[all.pc$trait=='RA',]

## identical

## examine the correlation between traits

pcor<-function(i){
  basis_trait<-names(ml)[i]
  message(basis_trait)
  bb_trait<-ml[[i]]
  title<-paste(basis_trait,bb_trait)
  basis.DT<-subset(DT,disease==basis_trait)[,.(id,lor,pwi)]
  basis.DT$label<-'basis'
  proj.DT<-subset(bb.DT,disease==bb_trait)[,.(id,lor,pwi)]
  proj.DT$label<-'proj'
  mDT<-rbind(basis.DT,proj.DT)
  mDT[,pwi_lor:=lor * pwi]
  mDT<-data.table::melt(mDT,id.vars=c('id','label'),measure.vars=c('lor','pwi_lor'))
  ## do lor
  mDT<-data.table::dcast(mDT,id~label+variable)
  par(mfrow=c(2,1))
  plot(mDT$basis_lor,mDT$proj_lor,main=basis_trait)
  plot(mDT$basis_pwi_lor,mDT$proj_pwi_lor,main=basis_trait)
}


mypdf = function(pdfname, mypattern = "MYTEMPPNG", ...) {
    fname = paste0(mypattern, "%05d.png")
    gpat = paste0(mypattern, ".*\\.png")
    takeout = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
    if (length(takeout) > 0)
        file.remove(takeout)
    pngname = file.path(tempdir(), fname)
    png(pngname, ...)
    return(list(pdfname = pdfname, mypattern = mypattern))
}
# copts are options to sent to convert
mydev.off = function(pdfname, mypattern, copts = "") {
    dev.off()
    gpat = paste0(mypattern, ".*\\.png")
    pngs = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
    mystr = paste(pngs, collapse = " ", sep = "")
    system(sprintf("convert %s -quality 100 %s %s", mystr, pdfname, copts))
}


res = mypdf("~/tmp/for_chris.pdf", res = 600, height = 7, width = 7, units = "in")
lapply(seq_along(ml),function(i){
  pcor(i)
})
mydev.off(pdfname = res$pdfname, mypattern = res$mypattern)
dev.off()



## taking RA as an example

ex<-subset(all.DT,disease %in% c('RA','bb:20002_1464:self_reported_rheumatoid_arthritis'))
setkey(ex,id)
ex<-ex[pwi]
ex[ex$disease=='bb:20002_1464:self_reported_rheumatoid_arthritis',disease:='bb_RA']
ex[,lp:=-log(p.val)]
## take the p.vals of the snps with top 10% pwi
n<-1
top.ex<-ex[ex$lp > quantile(ex$lp,prob=1-n/100),]

m<-melt(top.ex,id.vars=c('id','name','disease'),measure.vars='lor')
fin<-dcast(m,id+name~disease+variable)
fin<-subset(fin,!is.na(bb_RA_lor) & !is.na(RA_lor))

#library(cowplot)

ggplot(fin,aes(x=RA_lor,y=bb_RA_lor)) + geom_point() + theme_bw()
cor(fin$RA_lor,fin$bb_RA_lor)


## which id's don't seem to correlate
opp<-fin[which(fin$RA_lor>0 & fin$bb_RA_lor<0),]$id
## checked these and looks to be OK
