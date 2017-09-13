#library(data.table)
library(snpStats)
library(GenomicRanges)
library(Matrix)
library(corpcor)
library(mvtnorm)
library(magrittr)

## this script
## 1 creates simulations for aff.t1d
## 2 creates simulations for ill.t1d
## 3 uses these simulations to create confidence intervals for each metric which are integrated with
## those from simulate_AI_weightings_aff_ill.R and then saved.
setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")
source("./wakefield.R")
sim.size<-100
pi_1<-1e-4
#tmpfile<-"/scratch/ob219/as_basis/figure_data/analysis1_sims.RData"

# for ease of computation for simulation we will probably take forward gh_ss
# for computation without including quantitative traits
DT<-getGWASData()
## remove quantitative traits as well as meta.t1d
QT<-c('eosinophil','lymphocyte','myeloid','wbc')
DT<-DT[!DT$disease %in% c(QT,'meta.t1d'),]
## compute the different Metrics proposed
## beta == log(or) DT$lor already present
## Just use Z score
DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
## compute gamma - this can be computed in different ways
## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
##
DT[,gh_ss:=gamma_hat_ss(Z,total)]
DT[,gh_ss_pp:=gh_ss * pp]
# examine what happens if we compute gh using maf ?
DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
DT[,gh_maf_pp:=gh_maf * pp]
DT<-rbind(DT,createControl(DT))
## note that these files have previously been generated using ill.t1d - so data is already merged in
AT<-DT[,.(id,disease,ld.block,gh_ss,gh_maf,total,prop,Z,lor,maf,controls,cases)]

ST<-subset(AT,!disease %in% c('aff.t1d','ill.t1d'))

metrics<-c('lor','Z','gh_ss','gh_ss_pp','gh_maf','gh_maf_pp')
## this function generates a bunch of simulations for either aff.t1d or ill.t1d
#getSims<-function(sim.DT,,sim.size=sim.size){
getSims<-function(gwas=c('ill.t1d','aff.t1d'),sim_size=sim.size,under_null=FALSE){
  message(sprintf("Getting simulations for %s with %d simulations",gwas,sim.size))
  sim.DT<-subset(AT,disease==gwas)
  src_dir<-file.path(DATA_DIR,'as_basis/support/simulations/1KGenome_snpStats',gwas)
  snpStats.1kg.files<-list.files(path=src_dir,pattern="*.RData",full.names=TRUE)
	sim.ss.total<-unique(sim.DT$total)
	sim.ss.prop<-unique(sim.DT$prop)
	sim.by.chr<-lapply(snpStats.1kg.files,function(f){
	        message(sprintf("Processing %s",f))
	        # load in prepared snpMatrix objects - map contains annotations about SNPs gt contain genotypes
	        sm<-get(load(f))
	        ## precompute se_hat
	        sm$map$se_hat<-with(sm$map,lor/Z)
	        ## for some these will be NA estimate using above function
	        idx<-which(is.nan(sm$map$se_hat))
	        sm$map[idx,]$se_hat<-approxSE(sm$map[idx,]$maf,sim.ss.total)
					# need to add this so that GWASsim shrinks beta towards zero when the posterior for a SNP to have a non zero beta is low
					sm$map[,pp:=approx.bf.z2(z=Z,f=maf,N=sim.ss.total,s=sim.ss.prop,p=1e-4)]
	        GWASsim(sm,n.sims=sim.size,under.null=under_null)
	})
  all.chr.sim<-do.call('rbind',sim.by.chr) %>% melt(.,id.vars=c('id')) %>% merge(.,sim.DT,by.x='id',by.y='id')
  # compute posterior prob for inclusion
  all.chr.sim[,ppi:=approx.bf.z(value,maf,sim.ss.total,sim.ss.prop,pi_1),by=c('ld.block','variable')]
  ## reconstitute standard error
  all.chr.sim[,se:=lor/Z]
  all.chr.sim[is.nan(se),se:=approxSE(maf,sim.ss.total)]
  all.chr.sim[,Z:=value]
  all.chr.sim[,gh_ss:=gamma_hat_ss(Z,total)]
  all.chr.sim[,gh_ss_pp:=gh_ss * ppi]
  ## note here that Z has changed to the simulated version!
  all.chr.sim[,lor:=Z*se]
  # examine what happens if we compute gh using maf ?
  all.chr.sim[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(lor))]
  all.chr.sim[,gh_maf_pp:=gh_maf * ppi]
  setnames(all.chr.sim,'variable','simno')
  ## create a bunch of simulations for each of the different metrics
  ## create simulations
  sims<-lapply(metrics,function(m){
    message(sprintf("Processing %s",m))
    tmp<-melt(all.chr.sim,id.vars=c('id','simno'),measure.vars = m)
    tmp<-data.table::dcast(tmp,simno~id)
    tmp[,2:ncol(tmp),with=FALSE] %>% as.matrix(.)
  })
  names(sims)<-metrics
  return(sims)
}


# using 1kg we wish to simulate betas for ill.t1d study
#sim.DT<-subset(DT,disease=='ill.t1d')

## here we want to simulate aff for the projections and ill for the basis

ill.sim.file<-'/scratch/ob219/as_basis/figure_data/ill_sims_null.RData'
if(!file.exists(ill.sim.file)){
	ill.sims<-getSims('ill.t1d')
	save(ill.sims,file=ill.sim.file)
}else{
	load(ill.sim.file)
}

aff.sim.file<-'/scratch/ob219/as_basis/figure_data/aff_sims_null.RData'
if(!file.exists(aff.sim.file)){
	aff.sims<-getSims('aff.t1d',under_null=TRUE)
	save(aff.sims,file=aff.sim.file)
}else{
	load(aff.sim.file)
}

## compute outcomes of simulations for each of the different metrics
## load in actual data
(load("/scratch/ob219/as_basis/figure_data/analysis1.RData"))
## format so works with simulation dataset
## some simulated variants are zero (perhaps because lor = 0 remove these first - we assume that these are the same for all metrics


findMissingData<-function(sl){
	t<-sl[[1]]
	t<-rowSums(apply(t,1,is.na))
	names(t[t!=sim.size])
}

i.miss<-findMissingData(ill.sims)
a.miss<-findMissingData(aff.sims)
## take the intersection between two
keep<-intersect(i.miss,a.miss)

ill.sims<-lapply(ill.sims,function(s){
	s[,keep]
})

aff.sims<-lapply(aff.sims,function(s){
	s[,keep]
})

sim.bases<-lapply(bases,function(b){
	b[!rownames(b) %in% c('ill.t1d','aff.t1d'),colnames(b) %in% keep]
})

scaled.sims<-lapply(names(ill.sims),function(n){
	message(sprintf("Processing %s",n))
	sapply(1:sim.size,function(j){
     i.s<-ill.sims[[n]][j,]
     a.s<-aff.sims[[n]][j,]
     sbasis<-rbind(sim.bases[[n]],i.s)
	   pca<-prcomp(sbasis,center=TRUE,scale=FALSE)
	   ## project on another simulation
	   proj<-predict(pca,t(a.s))
	   vexp<-summary(pca)$importance[2,]
	   varWeightedEucledian(pca$x["i.s",],proj[1,],vexp)
	})
})
names(scaled.sims)<-names(ill.sims)

pw.dir<-'/home/ob219/scratch/as_basis/ai_ill_vs_aff_sims_null/'
pw.files <-list.files(path=pw.dir,pattern='*.RDS',full.names=TRUE)
pw.sims<-do.call('rbind',lapply(pw.files,readRDS))

scaled.sims$gh_ss_pw<-pw.sims[,'gh_ss_pw']
scaled.sims$gh_maf_pw<-pw.sims[,'gh_maf_pw']

unscaled.sims<-lapply(c('gh_ss_pp','gh_maf_pp'),function(sname){
	message(sprintf("Processing %s",sname))
        sapply(1:sim.size,function(j){
           i.s<-ill.sims[[sname]][j,]
	   	 		 # key difference here is that projection needs to be unscaled
	   			 bname<-sub('\\_p[pw]$','',sname)
           a.s<-aff.sims[[bname]][j,]
      		 sbasis<-rbind(sim.bases[[sname]],i.s)
      		 pca<-prcomp(sbasis,center=TRUE,scale=FALSE)
           ## project on another simulation
      		 proj<-predict(pca,t(a.s))
      		 vexp<-summary(pca)$importance[2,]
      		 varWeightedEucledian(pca$x["i.s",],proj[1,],vexp)
        })
})
names(unscaled.sims)<-paste('unscaled',c('gh_ss_pp','gh_maf_pp'),sep='_')

unscaled.sims$unscaled_gh_ss_pw<-pw.sims[,'gh_ss']
unscaled.sims$unscaled_gh_maf_pw<-pw.sims[,'gh_maf']



all.sims<-c(scaled.sims,unscaled.sims)
ci<-c(0.025,0.5,0.985)
empirical_confidence_intervals<-lapply(all.sims,quantile,probs=ci)
saveRDS(empirical_confidence_intervals,file='/home/ob219/scratch/as_basis/figure_data/analysis1_ill_vs_aff_ci_under_null.RDS')
