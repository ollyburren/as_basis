

## add argument processing

## this just simulates aff.t1d and then computes against the same basis

library(optparse)


option_list = list(
        make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file", metavar="character"),
        make_option(c("-s", "--sim_size"), type="integer", default=10,
                    help="Number of simulations to run. Note number run is doubled", metavar="character")
        )

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
if (is.null(args$out)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (output file)", call.=FALSE)
}

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

sim.size<-args$sim_size
pi_1<-1e-4
outfile<-args$out
metrics<-c('gh_ss_pw','gh_maf_pw')

# for ease of computation for simulation we will probably take forward gh_ss
# for computation
#without including quantitative traits
DT<-getGWASData()
## remove quantitative traits as well as meta.t1d
QT<-c('eosinophil','lymphocyte','myeloid','wbc')
DT<-DT[!DT$disease %in% c(QT,'meta.t1d'),]
## compute the different Metrics proposed
## beta == log(or) DT$lor already present
## Just use Z score
DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
# note here ill.t1d is a stub
#ST<-DT[!DT$disease %in% c('aff.t1d','ill.t1d'),]
## compute gamma - this can be computed in different ways
## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
##
DT[,gh_ss:=gamma_hat_ss(Z,total)]
# examine what happens if we compute gh using maf ?
DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
DT[,lp0:=log(1-approx.bf.z2(Z,maf,cases+controls,cases/(cases+controls),pi_1)),by=c('disease','ld.block')]
DT<-rbind(DT,createControl(DT),fill=TRUE)

## don't need all the columns
AT<-DT[,.(id,disease,lp0,ld.block,gh_ss,gh_maf,total,prop,Z,lor,maf,controls,cases)]

#ST<-subset(AT,!disease %in% c('aff.t1d','ill.t1d'))
#want to simulate this one only
ST<-subset(AT, disease!='aff.t1d')

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
	        GWASsim(sm,n.sims=sim_size,under.null=under_null)
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
}

aff.sims<-getSims('aff.t1d',under_null=TRUE)
#ill.sims<-getSims('ill.t1d')



## note that these files have previously been generated using ill.t1d - so data is already merged in

## a function to create paired DT for basis and projection
ST<-ST[ST$id %in% unique(aff.sims$id),]
aff.SF<-split(aff.sims,aff.sims$variable)

ST[ST$disease!='controls',q_i:=1-exp(sum(lp0)),by=id]
emp<-mean(ST[,list(mean_qi=mean(q_i)),by=id]$mean_qi)
## prior odds
po<-emp/(1-emp)
## note that emp us for h1 that beta != 0 therefore we need to take reciprocal as equation assumes pi_0 - see notes
po<-1/po
ST[,uABF:=po*(q_i/(1-q_i))]
BIG<-quantile(ST[is.finite(ST$uABF),]$uABF,prob=0.9999)
## BIG is still rather large but at least tractable
ST[is.infinite(uABF) | uABF > BIG, uABF:=BIG]
ST[,pwi:=computePWI(uABF,emp),by=c('disease','ld.block')]
ST[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]

getSims<-function(ST,sim){
  pwi<-unique(ST[,.(pwi),key='id'])
  setkey(sim,id)
  sim<-sim[pwi]
  sim[,c('gh_ss_pw','gh_maf_pw'):=list(gh_ss * pwi,gh_maf *pwi)]
  sim
}




## compute the basis for metrics we want
sims<-do.call('rbind',lapply(1:sim.size,function(i){
      message(sprintf("Processing %s",i))
      proj<-getSims(ST,aff.SF[[i]])
      do.call('c',lapply(metrics,function(m){
        tmp<-melt(ST,id.vars=c('id','disease'),measure.vars = m)
        basis.matrix<-data.table::dcast(tmp,id~disease)
        basis.matrix<-basis.matrix[,2:ncol(basis.matrix),with=FALSE] %>% as.matrix(.) %>% t(.)
        pca<-prcomp(basis.matrix,center=TRUE,scale=FALSE)
        ## now get the simulation convert to correct scale and project onto the basis
        res<-sapply(c(m,sub("_pw$","",m)),function(pm){
          tmp<-melt(proj,id.vars=c('id','disease'),measure.vars = pm)
          proj.v<-data.table::dcast(tmp,id~disease)
					proj.v<-proj.v[,2:ncol(proj.v),with=FALSE] %>% as.matrix(.) %>% t(.)
          proj<-predict(pca,proj.v)
					vexp<-summary(pca)$importance[2,]
					varWeightedEucledian(pca$x["ill.t1d",],proj[1,],vexp)
          })
      }))
}))


message(sprintf("Writing output to %s",outfile))
saveRDS(sims,file=outfile)
