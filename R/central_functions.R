library(data.table)
library(reshape2)
library(ggplot2)

DATA_DIR<-'/scratch/ob219/'
getGWASData<-function(){
	lf<-file.path(DATA_DIR,'as_basis/merged_data/17_traits.RData')
	message(sprintf("Loading trait data from: %s",lf))
	load(file.path(DATA_DIR,'as_basis/merged_data/17_traits.RData'))
	message(sprintf("Finished loading trait data from: %s",lf))
	## loads into final.t
	setkey(final.t,id)
	ss<-fread(file.path(DATA_DIR,'as_basis/gwas_stats/sample_counts.csv'))
	ss[,total:=cases+controls]
	ss[,prop:=signif(cases/(cases+controls),digits=2)]
	message("Merging study information")
	merge(final.t,ss,by.x='disease',by.y='disease')
}

createMatrix<-function(DT,var='lor'){
   message("Melting")
   DT<-melt(DT,id.vars=c('id','disease'),measure.vars = var)
   message("Recasting")
   ret<-data.table::dcast(DT,disease~id)
   diseases<-ret$disease
   ret<-as.matrix(ret[,2:ncol(ret),with=FALSE])
   rownames(ret)<-diseases
   return(ret)
   #message("Adding control")
   #fret<-rbind(ret,rep(0,ncol(ret)))
   #rownames(fret)<-c(diseases,'control')
   #return(fret)
}

# function computes PC for dataset without aff.t1d and then uses this to predict aff.t1d loadings
t1d.compare<-function(basis,center=TRUE,scale=FALSE){
	idx<-which(rownames(basis)=='aff.t1d')
	pca<-prcomp(basis[-idx,],center=center,scale=scale)
	proj<-predict(pca,t(basis[idx,]))
	list(pca=pca,proj=proj)
}

## compute weighted euc distance between PCA loadings
varWeightedEucledian<-function(simLoad,actLoad,vexp){
    sqrt((actLoad-simLoad)^2 %*% vexp)
}

## compares each trait variance weighted eucledian distance with cvar
compareVarWEuc<-function(pc,proj){
    	all.pc<-rbind(pc$x,proj)
   	 ## what is the variance explained
    	vexp<-summary(pc)$importance[2,]
    	apply(all.pc,1,varWeightedEucledian,all.pc[nrow(all.pc),],vexp)
}

gamma_hat_ss<-function(Z,total){
  rv<-numeric(length=length(Z))
  idx<-which(Z!=0)
  rv[idx]<-Z[idx] * (1/sqrt(2*total[idx]))
  return(rv)
}

ca<-function(n0,f){
    n0*(1-f)
}

cb<-function(n0,f){
    n0*f
}

cc<-function(n1,a,b,theta){
    (n1*a)/(a+(b*theta))
}

cd<-function(n1,a,b,theta){
    (n1*b)/(a+(b*theta))
}

gamma_hat_maf<-function(n0,n1,f,theta){
    n<-n0+n1
    a<-ca(n0,f)/n
    b<-cb(n0,f)/n
    c<-cc(n1,a,b,theta)/n
    d<-cd(n1,a,b,theta)/n
    recip.sm<-do.call('cbind',lapply(list(a,b,c,d),function(fi) 1/fi))
    return(sqrt(rowSums(recip.sm)))
}

createControl<-function(DT){
  ctrl<-subset(DT,!duplicated(DT$id))
  ctrl[,disease:='control']
  ctrl[,or:=1]
  ctrl[,p.val:=1]
  ctrl[,pp:=0]
  ctrl[,lor:=0]
  ctrl[,cases:=4000]
  ctrl[,controls:=4000]
  ctrl[,total:=cases+controls]
  ctrl[,prop:=0.5]
  ctrl[,Z:=0]
  ctrl[,gh_ss:=gamma_hat_ss(Z,total)]
  ctrl[,gh_ss_pp:=gh_ss * pp]
  # examine what happens if we compute gh using maf ?
  ctrl[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(or))]
  ctrl[,gh_maf_pp:=gh_maf * pp]
}

# This comes from R ggplot cookbook and is not my work
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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

GWASsim<-function(sm,n.sims=200,under.null=FALSE){
        map<-sm$map
        ld.idx<-split(1:nrow(map),map$LDBLOCK)
        by.chr<-lapply(seq_along(ld.idx),function(i){
                idx<-ld.idx[[i]]
                ## grab standard error estimates
                se_hat<-map[idx,]$se_hat
		# here we shrink our estimates for beta towards 0 when the posterior for a SNP to have a non zero beta is low
		# note because the standard error remains the same we still get Z scores that are not pp adjusted
		if(under.null){
			beta_hat<-rep(0,length(idx))
		}else{
			beta_hat<-map[idx,]$lor * map[idx,]$pp
		}
                if(length(idx)==1){
                        message("Only one sample from normal distro")
                        ## just sample from norm
                        ## var_beta is SE^2 but rnorm takes standard deviation so no need
                        ret<-t(rnorm(n.sims,mean=beta_hat,sd=se_hat))
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
                        ret<-mvs.perm(beta_hat,cov.beta,n=n.sims)
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


## this function computes a set of covariance matrices for use in GWASsim. 
## if we store these it is simple to simulate what happens with sample size by
## multiplying through covariance matrix by f(new.sample.size)/f(old.sample.size)

computeBetaCovariates<-function(sm){
	map<-sm$map
	ld.idx<-split(1:nrow(map),map$LDBLOCK)
        lapply(seq_along(ld.idx),function(i){
                idx<-ld.idx[[i]]
                ## grab standard error estimates
                se_hat<-map[idx,]$se_hat
                if(length(idx)==1){
                        message("Only one sample from normal distro")
			return(list(map=map[idx,],cov=se_hat))
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
			return(list(map=map[idx,],cov=cov.beta))
                }
        })
}

## this is a function to simulate shared GWAS 
GWASsimSampleSize<-function(i,scase=4000,sctrl=4000,ncase=4000,nctrl=4000,n.sims=200){
	## each one of these is chromosome
	ss.var.adj<-(1/(ncase+nctrl))/(1/(scase/sctrl))
	by.ld<-lapply(i,function(x){
		map<-x$map
		cov<-x$cov
		beta_hat<-map$lor * map$pp
		if(nrow(map)==1){
			ret<-t(rnorm(n.sims,mean=beta_hat,sd=cov * ss.var.adj^2))
		}else{
			ret<-mvs.perm(beta_hat,cov * ss.var.adj,n=n.sims)
		}
		(1/map$se_hat) * ret
	})
	## here we assume that the ld blocks are processed in chromosomal order
	ret<-data.table(do.call('rbind',by.ld))
	ret$id<-do.call('c',lapply(i,function(x) x$map$id))
	ret
}
