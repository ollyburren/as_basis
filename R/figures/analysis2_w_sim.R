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
tmpfile<-"/scratch/ob219/as_basis/figure_data/analysis2_sims.RData"

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
metrics<-c('lor','Z','gh_ss','gh_ss_pp','gh_maf','gh_maf_pp')
# using 1kg we wish to simulate betas for ill.t1d study
if(!file.exists(tmpfile)){
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
		## here we compute Z score simulations under the null
	        GWASsim(sm,under.null=TRUE)
	})
	all.chr.sim<-do.call('rbind',sim.by.chr) %>% melt(.,id.vars=c('id')) %>% merge(.,sim.DT,by.x='id',by.y='id')
	# compute posterior prob for inclusion
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
	metrics<-c('lor','Z','gh_ss','gh_ss_pp','gh_maf','gh_maf_pp')
	sims<-lapply(metrics,function(m){
		message(sprintf("Processing %s",m))
		tmp<-melt(all.chr.sim,id.vars=c('id','simno'),measure.vars = m)
		tmp<-data.table::dcast(tmp,simno~id)
		tmp[,2:ncol(tmp),with=FALSE] %>% as.matrix(.)
	})
	names(sims)<-metrics
	save(sims,file=tmpfile)
}else{
	load(tmpfile)
}

## compute outcomes of simulations for each of the different metrics
## load in actual data
(load("/scratch/ob219/as_basis/figure_data/analysis1.RData"))
## format so works with simulation dataset
## some simulated variants are zero (perhaps because lor = 0 remove these first - we assume that these are the same for all metrics
tsims<-sims[[1]]
tsims<-rowSums(apply(tsims,1,is.na))
sim.snps<-names(tsims[tsims!=200])
sims<-lapply(sims,function(s){
	s[,sim.snps]
})
sim.bases<-lapply(bases,function(b){
	b[!rownames(b) %in% c('aff.t1d'),colnames(b) %in% sim.snps]	
})
bases<-lapply(bases,function(b){
	b[,colnames(b) %in% sim.snps]	
})
# check that colnames are in the same order
#lapply(seq_along(sim.bases),function(i){
#	identical(colnames(sim.bases[[i]]),colnames(sims[[i]]))
#})
## simulate by computing pca with a simulation and then reprojecting another simulation
#scaled.sims<-lapply(seq_along(sim.bases),function(i){
#	message(sprintf("Processing %s",i))
#	sapply(1:(sim.size/2),function(j){
#	   sim<-sims[[i]][j,]
#	   p<-sims[[i]][j+(sim.size/2),]
#	   sbasis<-rbind(sim.bases[[i]],sim)
#	   pca<-prcomp(sbasis,center=TRUE,scale=FALSE)
#	   ## project on another simulation
#	   proj<-predict(pca,t(p))
#	   vexp<-summary(pca)$importance[2,]
#	   varWeightedEucledian(pca$x["sim",],proj[1,],vexp)
#	})
#})
cidx<-which(names(sim.bases) %in% c('gh_ss_pp','gh_maf_pp'))
unscaled.sims<-lapply(cidx,function(i){
	message(sprintf("Processing %s",i))
	## compute bases once
	pca<-prcomp(bases[[i]],center=TRUE,scale=FALSE)
        sapply(1:(sim.size/2),function(j){
	   # key difference here is that projection needs to be unscaled
	   bname<-sub('\\_pp$','',names(sim.bases)[i])
           p<-sims[[bname]][j+(sim.size/2),]
           ## project on another simulation
           proj<-predict(pca,t(p))
           vexp<-summary(pca)$importance[2,]
           varWeightedEucledian(pca$x["ill.t1d",],proj[1,],vexp)
        })
})
names(unscaled.sims)<-c('unscaled_gh_ss','unscaled_gh_maf')
ci<-c(0.025,0.5,0.985)
empirical_confidence_intervals_null<-lapply(unscaled.sims,quantile,probs=ci)
names(empirical_confidence_intervals_null)<-c('unscaled_gh_ss','unscaled_gh_maf')
## load in shared plot confidence intervals
all.sims.shared<-get(load("/scratch/ob219/as_basis/figure_data/analysis1_shared_simulation.RData"))
# Need to do unscaled input - then plot barplots using code from analysis1.R
empirical_confidence_intervals_shared<-lapply(all.sims.shared,quantile,probs=ci)
names(empirical_confidence_intervals_shared)<-c(metrics,c('unscaled_gh_ss','unscaled_gh_maf'))

## what about actual data ?
(load("/scratch/ob219/as_basis/figure_data/analysis1.RData"))
#scaled.comparisons<-lapply(bases,function(b){
#        comp<-t1d.compare(b)
#        with(comp,compareVarWEuc(pca,proj))
#})
#names(scaled.comparisons)<-metrics
# for gamma hat have a look at the effect of just using non ppi weightings
cidx<-which(names(bases) %in% c('gh_ss_pp','gh_maf_pp'))
unscaled.comparisons<-lapply(cidx,function(i){
        ## need to replace aff.t1d with unscaled version
        b<-bases[[i]]
        bname<-sub('\\_pp$','',names(bases)[i])
        t1d.gh<-DT[DT$disease=='aff.t1d',c(bname,'id'),with=FALSE]
        idx<-which(rownames(b)=='aff.t1d')
        ## get lors and make sure in the correct order
        b[idx,]<-t1d.gh[match(t1d.gh$id,colnames(b)),][[bname]]
        comp<-t1d.compare(b)
        with(comp,compareVarWEuc(pca,proj))
})

empirical_confidence_intervals_shared<-empirical_confidence_intervals_shared[which(names(empirical_confidence_intervals_shared) %in% c('unscaled_gh_ss','unscaled_gh_maf'))]

title<-list(expression("Unscaled"~hat(gamma[ss_ppi])),expression("Unscaled"~hat(gamma[maf_ppi])))
plots<-lapply(seq_along(unscaled.comparisons),function(i){
	mcin<-empirical_confidence_intervals_null[[i]]
	mshared<-empirical_confidence_intervals_shared[[i]]
        m<-names(unscaled.comparisons)[i]
        dat<-unscaled.comparisons[[i]]
        dat<-dat[dat!=0]
        dat<-data.frame(disease=names(dat),distance=dat)
        ## hilight ill.t1d and control sets
        dat$hilight<-logical(length=nrow(dat))
        dat[dat$disease %in% c('ill.t1d','control'),]$hilight<-TRUE
        dat$disease<-factor(dat$disease,levels=dat[order(dat$distance),]$disease)
        ggplot(dat,aes(x=disease,y=distance,fill=hilight,color=hilight)) + geom_bar(stat="identity") + theme_bw() + xlab("Trait") + ylab("Distance") + ggtitle(title[[i]]) + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + scale_fill_manual(name="Bar",values=c('white','firebrick'),guide=FALSE) + scale_color_manual(name="Bar",values=c('black','black'),guide=FALSE) + 
geom_rect(xmin=0,xmax=length(levels(dat$disease)) + 0.5,ymin=mcin[1],ymax=mcin[3],alpha=0.1,color=NA,fill='grey90') + 
geom_hline(yintercept=mcin[2],color='firebrick') + 
geom_rect(xmin=0,xmax=length(levels(dat$disease)) + 0.5,ymin=mshared[1],ymax=mshared[3],alpha=0.1,color=NA,fill='steelblue') + 
geom_hline(yintercept=mshared[2],color='blue')
})

pdf("/home/ob219/git/as_basis/pdf/analysis2_w_sim.pdf")
multiplot(plotlist=plots,cols=1)
dev.off()
