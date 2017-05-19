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
   message("Adding control")
   fret<-rbind(ret,rep(0,ncol(ret)))
   rownames(fret)<-c(diseases,'control')
   return(fret)
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
