## this software creates a figure that compares the eucledian distances for different scaling metrics
#run this from here
setwd("/home/ob219/git/as_basis/R/")
## 1. import data and compute different scalings
source("./central_functions.R")
tmpfile<-"/scratch/ob219/as_basis/figure_data/analysis1_w_QT.RData"
if(file.exists(tmpfile)){
        load(tmpfile)
}else{
	DT<-getGWASData()
	## remove quantitative traits as well as meta.t1d
	DT<-DT[DT$disease != 'meta.t1d',]
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
	#computing gamma hat using MAF does not make sense for QT
	metrics<-c('lor','Z','gh_ss','gh_ss_pp')
	bases<-lapply(metrics,function(m){
        	createMatrix(DT,var=m)
	})
	names(bases)<-metrics
	save(bases,file=tmpfile)
}


## plots for scaled comparisons
scaled.comparisons<-lapply(bases,function(b){
	comp<-t1d.compare(b)
	with(comp,compareVarWEuc(pca,proj))
})
names(scaled.comparisons)<-metrics
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

all.comparisons<-c(scaled.comparisons,unscaled.comparisons)
title<-list(expression(hat(beta)),expression(Z),expression(hat(gamma[ss])),expression(hat(gamma[ss_ppi])),expression("Non PPi Input Projection"~hat(gamma[ss_ppi])))
plots<-lapply(seq_along(all.comparisons),function(i){
        m<-names(all.comparisons)[i]
        dat<-all.comparisons[[i]]
        dat<-dat[dat!=0]
        dat<-data.frame(disease=names(dat),distance=dat)
	## hilight ill.t1d and control sets
	dat$hilight<-logical(length=nrow(dat))
	dat[dat$disease %in% c('ill.t1d','control'),]$hilight<-TRUE
	dat$disease<-factor(dat$disease,levels=dat[order(dat$distance),]$disease)
	ggplot(dat,aes(x=disease,y=distance,fill=hilight,color=hilight)) + geom_bar(stat="identity") + theme_bw() + xlab("Trait") + ylab("Distance") + ggtitle(title[[i]]) + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + scale_fill_manual(name="Bar",values=c('white','firebrick'),guide=FALSE) + scale_color_manual(name="Bar",values=c('black','black'),guide=FALSE)
})
pdf("/home/ob219/git/as_basis/pdf/analysis1_w_QT.pdf")	
multiplot(plotlist=plots,cols=3)
dev.off()
