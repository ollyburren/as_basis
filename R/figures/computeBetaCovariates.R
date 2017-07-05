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
sim.size<-200
pi_1<-1e-4
tmpfile<-"/scratch/ob219/as_basis/figure_data/cov_matrix_ill.t1d.RData"

# for ease of computation for simulation we will probably take forward gh_ss
# for computation without including quantitative traits
DT<-getGWASData()
## remove quantitative traits as well as meta.t1d
QT<-c('eosinophil','lymphocyte','myeloid','wbc')
DT<-DT[DT$disease == 'ill.t1d',]
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
## note that these files have previously been generated using ill.t1d - so data is already merged in
snpStats.1kg.files<-list.files(path=file.path(DATA_DIR,'as_basis/support/simulations/1KGenome_snpStats'),pattern="*.RData",full.names=TRUE)
## here we compute beta ~ MVN(0,Sigma) - this simulates where there is no overlap (i.e. no signal) not sure that this is completely what we want to do.
sim.DT<-subset(DT,disease=='ill.t1d')
sim.ss.total<-unique(sim.DT$total)
sim.ss.prop<-unique(sim.DT$prop)
# using 1kg we wish to simulate betas for ill.t1d study
if(!file.exists(tmpfile)){
	cov.by.chr<-lapply(snpStats.1kg.files,function(f){
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
		## here we compute Z score simulations under the null
		computeBetaCovariates(sm)
	})
	## here save sim.DT and cov matrix
	save(list=c('cov.by.chr','sim.DT'),file=tmpfile)
}else{
	message(sprinf("File %s already exists. Remove to refresh !"),tmpfile)
}
