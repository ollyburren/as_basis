setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")
source("./wakefield.R")
library(magrittr)
sim.size<-100
pi_1<-1e-4
tmpfile<-"/home/ob219/scratch/as_basis/merged_data/AI_raw.RDS"

# for ease of computation for simulation we will probably take forward gh_ss
# for computation without including quantitative traits

if(!file.exists(tmpfile)){
    DT<-getGWASData()
    ## remove quat traits
    QT<-c('eosinophil','lymphocyte','myeloid','wbc')
    RMT<-c(QT,'aff.t1d','ill.t1d')
    jia<-c('jia_cc','jia_EO','jia_ERA','JIA_nosys','jia_PO','jia_PsA','jia_RFneg','jia_RFpos','jia_sys')
    DT<-DT[!DT$disease %in% RMT,]

    DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
    ## compute gamma - this can be computed in different ways
    ## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
    ##
    DT[,gh_ss:=gamma_hat_ss(Z,total)]
    DT[,gh_ss_pp:=gh_ss * pp]
    # examine what happens if we compute gh using maf ?
    DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
    DT[,gh_maf_pp:=gh_maf * pp]
    jia.DT<-DT[DT$disease %in% jia,]
    DT<-DT[!DT$disease %in% jia,]
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

jia.DT<-getGWASData()
## remove quat traits
jia<-c('jia_cc','jia_EO','jia_ERA','JIA_nosys','jia_PO','jia_PsA','jia_RFneg','jia_RFpos','jia_sys')
jia.DT<-jia.DT[jia.DT$disease %in% jia,]

jia.DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
## compute gamma - this can be computed in different ways
## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
##
jia.DT[,gh_ss:=gamma_hat_ss(Z,total)]
jia.DT[,gh_ss_pp:=gh_ss * pp]
# examine what happens if we compute gh using maf ?
jia.DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
jia.DT[,gh_maf_pp:=gh_maf * pp]
setkey(jia.DT,'id')
## need to add weightings
pwi<-unique(DT[,.(id,pwi),key='id'])
jia.DT<-jia.DT[pwi]
jia.DT[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]


## metric/weighting to use
metric<-'gh_maf_pw'
## trait to project
ptrait<-'jia_cc'

asb<-createBasis(DT,metric)

## checked and the order is the same between basis and proj
proj<-subset(jia.DT,disease==ptrait)[[metric]] %>% as.matrix(.) %>% t(.)

pred<-predict(asb,proj)

p <- autoplot(asb, label = TRUE, label.size = 3)
p + geom_point(x=pred[1,1],y=pred[1,2],color='red')

## work out how to draw scree plots to compare at each PC
## also work out how to use euclidian stuff to compare.
