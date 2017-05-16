library(data.table)
library(snpStats)
library(GenomicRanges)
library(Matrix)
library(corpcor)
library(mvtnorm)
library(reshape2)
library(ggplot2)


## FUNCTION

## if we have OR 1 then Z is 0 we can approximate SE using the following
approxSE<-function(f,N){
        (1/sqrt(f*(1-f))) * (1/sqrt(2*N))
}

## Code to sample multivariate norm
mvs.perm<-function(m,sigma,n=1000){
        if(!is.matrix(sigma))
                stop("sigma parameter is not a matrix")
        if(!is.positive.definite(sigma,,method="chol"))
                stop("sigma is not positive definite")
        rd<-rmvnorm(n,mean=m,sigma=sigma,method="chol")
        t(rd)
}


##  Code to compute sigma - genotype covariance matrix
mvs.sigma.r<-function(r){
        diag(r)<-1
        if(any(is.na(r)))
                message(sprintf("Found %s where r was na",sum(is.na(r))))
                r[is.na(r)]<-0
        return(as(make.positive.definite(r),"Matrix"))
}

# simulate a GWAS based on log(OR) from an actual GWAS
GWASsim<-function(sm,n.sims=200){
	map<-sm$map
        ld.idx<-split(1:nrow(map),map$LDBLOCK)
        by.chr<-lapply(seq_along(ld.idx),function(i){
                idx<-ld.idx[[i]]
                ## grab standard error estimates
                se_hat<-map[idx,]$se_hat
		beta_hat<-map[idx,]$lor
                if(length(idx)==1){
                        message("Only one sample from normal distro")
                        ## just sample from norm
                        ## var_beta is SE^2 but rnorm takes standard deviation so no need
                        ret<-t(rnorm(200,mean=beta_hat,sd=se_hat))
                }else{
                        #multivariate
                        gt<-sm$gt[,idx]
                        colnames(gt)<-map[idx,]$POS
			# compute R statistic
                        r<-ld(gt,gt,stats="R")
			# compute closest pos-def covariance matrix
                        r<-as.matrix(mvs.sigma.r(Matrix(r)))
                        ## for beta the covariance matrix is estimates by sigma x SE * SE^T
                        cov_se<-tcrossprod(se_hat)
                        cov.beta<-cov_se * r
                        ## simulate beta
                        ret<-mvs.perm(beta_hat,cov.beta,n=200)
                }
		# compute Z score
                (1/se_hat) * ret
        })
        ret<-data.table(do.call('rbind',by.chr))
        #ret$id<-sm$map$id
	# perhaps safer
	ret$id<-do.call('c',split(map$id,map$LDBLOCK))
        ret	
}

## turn a long thin data table into a shorter fatter one based on a particular variable
createMatrix<-function(DT,var='lor'){
   DT<-melt(DT,id.vars=c('id','disease'),measure.vars = var)
   ret<-data.table::dcast(DT,disease~id)
   diseases<-ret$disease
   ret<-as.data.frame(ret[,2:ncol(ret),with=FALSE])
   rownames(ret)<-diseases
   fret<-rbind(ret,rep(0,ncol(ret)))
   rownames(fret)<-c(diseases,'control')
   fret
}

## compute weighted euc distance between PCA loadings
nEucledian<-function(simLoad,actLoad,vexp){
    sqrt((actLoad-simLoad)^2 %*% vexp)
}

# compute weighted Euc distance between proj and simulated PC only
computeVarWEuc<-function(pc,proj.pc){
    all.pc<-rbind(pc$x,proj.pc)
    ## what is the variance explained
    vexp<-summary(pc)$importance[2,]
    idx<-nrow(all.pc)
    apply(all.pc,1,nEucledian,all.pc[idx,],vexp)["sim"]
}

# compute weighted Euc distance between proj and actual basis traits
compareVarWEuc<-function(pc,proj.pc,cvar){
    all.pc<-rbind(pc$x,proj.pc)
    ## what is the variance explained
    vexp<-summary(pc)$importance[2,]
    idx<-which(rownames(all.pc)==cvar)
    apply(all.pc,1,nEucledian,all.pc[idx,],vexp)
}

DATA_DIR<-'/scratch/ob219/'
## load processed summary statistics for all 17 traits
(load(file.path(DATA_DIR,'as_basis/merged_data/17_traits.RData')))
## loads into final.t 
setkey(final.t,id)
target.gwas<-'ill.t1d'
## estimate z score from p.value
final.t$Z<-qnorm(0.5 * final.t$p.val, lower.tail = FALSE) * sign(final.t$lor)
of<-list.files(path=file.path(DATA_DIR,'as_basis/support/simulations/1KGenome_snpStats'),pattern="*.RData",full.names=TRUE)
#load in sample size information about each study
ss<-fread(file.path(DATA_DIR,'as_basis/gwas_stats/sample_counts.csv'))
ss[,total:=cases+controls]
ss[,prop:=signif(cases/(cases+controls),digits=2)]
target.gwas.total.sample.size<-subset(ss,disease==target.gwas)$total
target.gwas.prop.case.control<-subset(ss,disease==target.gwas)$prop


## for each chromosome and LD block create 200 simulations
sim.by.chr<-lapply(of,function(f){
	message(sprintf("Processing %s",f))
	# load in prepared snpMatrix objects - map contains annotations about SNPs gt contain genotypes
	sm<-get(load(f))
	## precompute se_hat
	sm$map$se_hat<-with(sm$map,lor/Z)
	## for some these will be NA estimate using above function
	idx<-which(is.nan(sm$map$se_hat))
	sm$map[idx,]$se_hat<-approxSE(sm$map[idx,]$maf,ss[ss$disease==target.gwas,]$total)
	GWASsim(sm)
})
# OK so have a matrix of Z scores
all.chr.sim<-do.call('rbind',sim.by.chr)
## note that some SNPs are missing from the simulations. Not sure why but unlikely 
## to affect things too much.
all.ids<-do.call('c',lapply(sim.by.chr,'[[','id'))
tmp<-melt(all.chr.sim,id.vars=c('id'))
## compute pAdjZ (the posterior prob adjusted Z score ppi * Z * 1/sqrt(N)
to.comp.ppi<-final.t[final.t$disease==target.gwas,c('ld.block','maf','id'),with=FALSE]
tmpy<-merge(tmp,to.comp.ppi,by.x='id',by.y='id')
## try and d  things the data.table way
## compute ppi
source(file.path(DATA_DIR,'as_basis/tmp/wakefield.R'))
tmpy[,ppi:=approx.bf.z(value,maf,target.gwas.total.sample.size,target.gwas.prop.case.control,1e-4),by=c('ld.block','variable')]
## constant variance due to sample size
var_sam<-1/sqrt(2*(target.gwas.total.sample.size))
tmpy$nv<-tmpy$value * var_sam * tmpy$ppi
tmpy$value<-tmpy$nv
tmpy<-tmpy[,1:3,with=FALSE]
tmpy.sim<-dcast(tmpy,variable~id)
tmpy.sim<-tmpy.sim[,2:ncol(tmpy.sim)]
## columns here are different to those in basis need to reorder
ss<-ss[,c('disease','total'),with=FALSE]
final.t.ss<-merge(final.t,ss,by.x='disease',by.y='disease')
#multiply by two to get allele counts
final.t.ss$Zadj<-final.t.ss$Z * sqrt(1/(2*final.t.ss$total))
final.t.ss$pZadj<-final.t.ss$Zadj * final.t.ss$pp 
##there are some variants missing - for this analysis it should not matter too much as long as it's not an error.
## instead we need to make sure that these agree with the actual scores we are using for the pca basis.
filt<-subset(final.t.ss,id %in% all.chr.sim$id & disease %in% c('asthma','CD','CEL','eosinophil','JIA_nosys','lymphocyte','MS','myeloid','PBC','PSC','RA','SLE','UC','wbc'))
## make sure that the same variants are in simulations as in basis dataset
setdiff(all.chr.sim$id,filt$id)
## appears OK

## create the matrix
basis<-createMatrix(filt,var='pZadj')
setcolorder(tmpy.sim,colnames(basis))
basis<-as.matrix(basis)
tmpy.sim<-as.matrix(tmpy.sim)
## next we add a row from sim and compute pca then project from proj
## will take approx an hour to compute serially
results<-lapply(1:100,function(i){
	message(sprintf("Processing %d",i))
	mat<-rbind(basis,tmpy.sim[i,])
	rownames(mat)[length(rownames(mat))]<-'sim'
	pca<-prcomp(mat,center=TRUE,scale=FALSE)
	proj.pc<-predict(pca,t(tmpy.sim[i+100,]))
	computeVarWEuc(pca,proj.pc)
})

## lets look at the variance
r<-do.call('c',results)
v<-numeric(length=length(r)-1)
for(i in 2:length(r)) v[i-1]<-var(r[1:i])
ggplot(data.table(index=1:length(v),variance=v),aes(x=index,y=variance)) + geom_point()
## not sure that it constitutes convergence but let's go with this ballpark for time being
## create true basis and conduct analysis

## do this again to get correct data
act<-subset(final.t, id %in% all.ids & disease %in% c('asthma','CD','CEL','eosinophil','JIA_nosys','lymphocyte','MS','myeloid','PBC','PSC','RA','SLE','UC','wbc','ill.t1d','aff.t1d') )
## estimate z score from p.value
act$Z<-qnorm(0.5 * act$p.val, lower.tail = FALSE) * sign(act$lor)
## need to recalculate Zadj for sample size
act<-merge(act,ss,by.x='disease',by.y='disease')
act$Zadj<-act$Z * 1/sqrt(2*act$total)
act$pZadj<-act$Zadj * act$pp 
act.basis<-createMatrix(act,var='pZadj')
## remove aff.t1d that we want to predict
idx<-which(rownames(act.basis)=='aff.t1d')
pca.act<-prcomp(act.basis[-idx,],center=TRUE,scale=FALSE)
## project aff.t1d onto basis
proj.aff.t1d<-predict(pca.act,act.basis[idx,])
sc.euc<-compareVarWEuc(pca.act,proj.aff.t1d,'aff.t1d')
DT<-data.table(disease=names(sc.euc),scaled=sc.euc)
DT$disease<-factor(DT$disease,levels = DT[order(DT$scaled),]$disease)
ggplot(DT,aes(x=disease,y=scaled)) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle = -90, hjust = 0))
