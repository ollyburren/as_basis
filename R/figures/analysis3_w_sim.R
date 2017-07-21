library(optparse)


option_list = list(
  	make_option(c("-c", "--cases"), type="numeric", default=4000,
              help="Number of cases to simulate", metavar="numeric"),
  	make_option(c("-t", "--controls"), type="numeric", default=4000,
              help="Number of controls to simulate", metavar="numeric"),
	make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$out)){
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
sim.size<-100
no.ctrl<-4000
no.case<-4000
## when deploying these will be parameterised so we can see what happens with different case control configurations
sim.case<-opt$cases
sim.ctrl<-opt$controls
pi_1<-1e-4
cov.file<-"/scratch/ob219/as_basis/figure_data/cov_matrix_ill.t1d.RData"
load(cov.file)
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
snpStats.1kg.files<-list.files(path=file.path(DATA_DIR,'as_basis/support/simulations/1KGenome_snpStats'),pattern="*.RData",full.names=TRUE)
## here we compute beta ~ MVN(0,Sigma) - this simulates where there is no overlap (i.e. no signal) not sure that this is completely what we want to do.
sim.DT<-subset(DT,disease=='ill.t1d')
sim.ss.total<-unique(sim.DT$total)
sim.ss.prop<-unique(sim.DT$prop)
## only want to simulate what happens for these scalings as these appear to be the most promising
metrics<-c('gh_ss_pp','gh_maf_pp')
## load in beta that we have already computed

processSims<-function(SIMS,metrics){
	all.chr.sim<-do.call('rbind',SIMS) %>% melt(.,id.vars=c('id')) %>% merge(.,sim.DT,by.x='id',by.y='id')
	all.chr.sim[,ppi:=approx.bf.z(value,maf,sim.ss.total,sim.ss.prop,pi_1),by=c('ld.block','variable')]
	## reconstitute standard error
	all.chr.sim[,se:=lor/Z]
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
	tsims<-sims[[1]]
	tsims<-rowSums(apply(tsims,1,is.na))
	sim.snps<-names(tsims[tsims!=200])
	lapply(sims,function(s){
        	s[,sim.snps]
	})
}

## use these simulations to create the basis - therefore they have a constant sample size (i.e. the adj value =1)
basis.sims<-lapply(cov.by.chr,GWASsimSampleSize,scase=4000,sctrl=4000,nctrl=4000,ncase=4000,n.sims=sim.size)
basis.sims<-processSims(basis.sims,c('gh_ss_pp','gh_maf_pp'))


sims<-lapply(cov.by.chr,GWASsimSampleSize,scase=sim.case,sctrl=sim.ctrl,nctrl=no.ctrl,ncase=no.case,n.sims=sim.size)
sims<-processSims(sims,c('gh_ss','gh_maf'))
## compute outcomes of simulations for each of the different metrics
## load in actual data goes into 'bases'
## format so works with simulation dataset
## some simulated variants are zero (perhaps because lor = 0 remove these first - we assume that these are the same for all metrics
(load("/scratch/ob219/as_basis/figure_data/analysis1.RData"))
bases<-lapply(bases,function(b){
	b[!rownames(b) %in% c('aff.t1d','ill.t1d'),colnames(b) %in% colnames(sims[[1]])]
})
cidx<-which(names(basis.sims) %in% c('gh_ss_pp','gh_maf_pp'))
unscaled.sims<-lapply(names(basis.sims),function(n){
        message(sprintf("Processing %s",n))
        sapply(1:sim.size,function(j){
	   for.basis<-basis.sims[[n]][j,]
           # key difference here is that projection needs to be unscaled
           bname<-sub('\\_pp$','',n)
           for.proj<-sims[[bname]][j,]
	   ## get precomputed bases
	   cp.basis<-bases[[n]]
	   ## some of these might be NA need to remove
	   idx<-which(!is.na(for.basis))
	   for.basis<-for.basis[idx]
	   cp.basis<-cp.basis[,idx]
           sbasis<-rbind(cp.basis,for.basis)
           pca<-prcomp(sbasis,center=TRUE,scale=FALSE)
           ## project on another simulation
           proj<-predict(pca,t(for.proj[idx]))
           vexp<-summary(pca)$importance[2,]
           varWeightedEucledian(pca$x["for.basis",],proj[1,],vexp)
        })
})
names(unscaled.sims)<-c('unscaled_gh_ss','unscaled_gh_maf')
ci<-c(0.025,0.5,0.985)
empirical_confidence_intervals_sim<-lapply(unscaled.sims,quantile,probs=ci)
names(empirical_confidence_intervals_sim)<-c('unscaled_gh_ss','unscaled_gh_maf')
## load in shared plot confidence intervals
obj<-list(case.no=opt$cases,ctrl.no=opt$controls,sims=unscaled.sims,ci=empirical_confidence_intervals_sim)
save(obj,file=opt$out)
message(sprintf("Saved results for %d cases and %d controls to %s",obj$case.no,obj$ctrl.no,opt$out))
