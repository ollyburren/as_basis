## this software creates a figure that compares the eucledian distances for different scaling metrics
#run this from here
setwd("/home/ob219/git/as_basis/R/")
## 1. import data and compute different scalings
source("./central_functions.R")
DT<-getGWASData()
## remove meta.t1d
DT<-DT[DT$disease!='meta.t1d',]
## compute the different Metrics proposed
## beta == log(or) DT$lor already present
## Just use Z score
DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
## compute gamma - this can be computed in different ways
## by using number of cases and controls if we were to include just case control in the basis - for time being just stick with this
## 
DT[,gh_ss:=gamma_hat_ss(Z,total)]
DT[,gh_ss_pp:=gh_ss * pp]
metrics<-c('lor','Z','gh_ss','gh_ss_pp')
bases<-lapply(metrics,function(m){
        createMatrix(DT,var=m)
})
names(bases)<-metrics


## plots for scaled comparisons
scaled.comparisons<-lapply(bases,function(b){
	comp<-t1d.compare(b)
	with(comp,compareVarWEuc(pca,proj))
})
names(scaled.comparisons)<-metrics
# for gamma hat have a look at the effect of just using non ppi weightings
cidx<-which(names(bases) %in% c('gh_ss_pp'))
t1d.gh<-DT[DT$disease=='aff.t1d',c('gh_ss','id'),with=FALSE]
unscaled.comparisons<-lapply(bases[cidx],function(b){
	## need to replace aff.t1d with unscaled version
	idx<-which(rownames(b)=='aff.t1d')
	## get lors and make sure in the correct order
	b[idx,]<-t1d.gh[match(t1d.gh$id,colnames(b)),]$gh
	comp<-t1d.compare(b)
        with(comp,compareVarWEuc(pca,proj))
})

all.comparisons<-c(scaled.comparisons,unscaled.comparisons)
title<-list(expression(hat(beta)),expression(Z),expression(hat(gamma)),expression(hat(gamma[ppi])),expression("Non PPi Input Projection"~hat(gamma[ppi])))
plots<-lapply(seq_along(all.comparisons),function(i){
        m<-names(all.comparisons)[i]
        dat<-all.comparisons[[i]]
        dat<-dat[dat!=0]
        dat<-data.frame(disease=names(dat),distance=dat)
	dat$disease<-factor(dat$disease,levels=dat[order(dat$distance),]$disease)
	ggplot(dat,aes(x=disease,y=distance)) + geom_bar(stat="identity") + theme_bw() + xlab("Trait") + ylab("Distance") + ggtitle(title[[i]]) + theme(axis.text.x=element_text(angle = -90, hjust = 0))
})
pdf("/home/ob219/git/as_basis/pdf/analysis1.pdf")	
multiplot(plotlist=plots,cols=2)
dev.off()
