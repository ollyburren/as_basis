library(data.table)
library(snpStats)
library(GenomicRanges)
library(Matrix)
library(corpcor)
library(mvtnorm)
library(magrittr)

setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")
source("./wakefield.R")
sim.size<-10
pi_1<-1e-4
tmpfile<-"/scratch/ob219/as_basis/figure_data/test_analysis1_sims.RData"
metrics<-c('gh_ss_pw','gh_maf_pw')

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
# note here ill.t1d is a stub
ST<-DT[!DT$disease %in% c('aff.t1d','ill.t1d'),]
## compute gamma - this can be computed in different ways
## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
##
ST[,gh_ss:=gamma_hat_ss(Z,total)]
# examine what happens if we compute gh using maf ?
ST[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
ST[,lp0:=log(1-approx.bf.z2(Z,maf,cases+controls,cases/(cases+controls),pi_1)),by=c('disease','ld.block')]
ST<-rbind(ST,createControl(ST),fill=TRUE)

## don't need all the columns
ST<-ST[,c('id','disease','lp0','ld.block','gh_ss','gh_maf'),with=FALSE]

## note that these files have previously been generated using ill.t1d - so data is already merged in
snpStats.1kg.files<-list.files(path=file.path(DATA_DIR,'as_basis/support/simulations/1KGenome_snpStats'),pattern="*.RData",full.names=TRUE)
# using 1kg we wish to simulate betas for ill.t1d study
sim.DT<-subset(DT,disease=='ill.t1d')
sim.ss.total<-unique(sim.DT$total)
sim.ss.prop<-unique(sim.DT$prop)
#if(!file.exists(tmpfile)){
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
	        GWASsim(sm,n.sims=sim.size)
	})
	all.chr.sim<-do.call('rbind',sim.by.chr) %>% melt(.,id.vars=c('id')) %>% merge(.,sim.DT,by.x='id',by.y='id')
	all.chr.sim[,se:=lor/Z]
	## some of these will be zero so Approximate
	all.chr.sim[is.nan(se),se:=approxSE(maf,sim.ss.total)]
	all.chr.sim[,Z:=value]
	all.chr.sim[,gh_ss:=gamma_hat_ss(Z,total)]
	## note here that Z has changed to the simulated version!
	all.chr.sim[,lor:=Z*se]
	# examine what happens if we compute gh using maf ?
	all.chr.sim[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(lor))]
	all.chr.sim[,lp0:=log(1-approx.bf.z2(Z,maf,cases+controls,cases/(cases+controls),pi_1)),by=c('variable','ld.block')]
	## after simulation some of SNPs are removed or cannot be simulated these need to be removed from the basis
	## do some tidying up as don't need all columns and just make things slow
	all.chr.sim<-all.chr.sim[,c('id','disease','lp0','ld.block','gh_ss','gh_maf','variable'),with=FALSE]

	## a function to create paired DT for basis and projection
	ST<-ST[ST$id %in% unique(all.chr.sim$id),]
	SF<-split(all.chr.sim,all.chr.sim$variable)
	getMatchedSims<-function(ST,bsim,psim){
		sits<-rbind(ST,bsim,fill=TRUE)
		## compute lpo for simulation
		sits[,q_i:=1-exp(sum(lp0)),by=id]
		## our empirical prior for these diseases is then
		emp<-mean(sits[,list(mean_qi=mean(q_i)),by=id]$mean_qi)
		po<-emp/(1-emp)
		## prior odds
		po<-emp/(1-emp)
		## note that emp us for h1 that beta != 0 therefore we need to take reciprocal as equation assumes pi_0 - see notes
		po<-1/po
		sits[,uABF:=po*(q_i/(1-q_i))]
		BIG<-quantile(sits[is.finite(sits$uABF),]$uABF,prob=0.9999)
		## BIG is still rather large but at least tractable
		sits[is.infinite(uABF) | uABF > BIG, uABF:=BIG]
		sits[,pwi:=computePWI(uABF,emp),by=c('disease','ld.block')]
		pwi<-unique(sits[,.(pwi),key='id'])
		setkey(psim,id)
		psim<-psim[pwi]
		psim[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]
		sits[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]
		list(basis=sits,proj=psim)
	}

	matched.sims<-do.call('rbind',lapply(1:(sim.size/2),function(i){
		message(sprintf("Processing %s",i))
		ms<-getMatchedSims(ST,SF[[i]],SF[[(sim.size/2)+i]])
		do.call('c',lapply(c('gh_ss_pw','gh_maf_pw'),function(m){
				bmt<-ms[[1]]
				proj<-ms[[2]]
				tmp<-melt(bmt,id.vars=c('id','disease'),measure.vars = m)
				basis.matrix<-data.table::dcast(tmp,id~disease)
				basis.matrix<-basis.matrix[,2:ncol(basis.matrix),with=FALSE] %>% as.matrix(.) %>% t(.)
				pca<-prcomp(basis.matrix,center=TRUE,scale=FALSE)
				res<-sapply(c(m,sub("_pw$","",m)),function(pm){
					tmp<-melt(proj,id.vars=c('id','disease'),measure.vars = pm)
					proj.v<-data.table::dcast(tmp,id~disease)
					proj.v<-proj.v[,2:ncol(proj.v),with=FALSE] %>% as.matrix(.) %>% t(.)
					## i checked and column id's are both the same in projection and basis
					## project on another simulation
					proj<-predict(pca,proj.v)
					vexp<-summary(pca)$importance[2,]
					varWeightedEucledian(pca$x["ill.t1d",],proj[1,],vexp)
				})
		}))
}))
