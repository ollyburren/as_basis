## this software creates a figure that compares the eucledian distances for different scaling metrics
#run this from here
setwd("/home/ob219/git/as_basis/R/")
## 1. import data and compute different scalings
source("./central_functions.R")
source("./wakefield.R")
pi_1<-1e-4
metrics<-c('lor','Z','gh_ss','gh_ss_pp','gh_ss_pw','gh_maf','gh_maf_pp','gh_maf_pw')
tmpfile<-"/scratch/ob219/as_basis/figure_data/analysis1.RData"
if(file.exists(tmpfile)){
	load(tmpfile)
}else{
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
	DT[,gh_maf:=lor/gamma_hat_maf(controls,cases,maf,exp(lor))]
	DT[,gh_maf_pp:=gh_maf * pp]
	DT<-rbind(DT,createControl(DT))
	##
	DT.no.affy<-DT[DT$disease != 'aff.t1d']
	DT.no.affy[,lp0:=log(1-approx.bf.z2(Z,maf,cases+controls,cases/(cases+controls),pi_1)),by=c('disease','ld.block')]
	## compute q_i across all diseases for i-th SNP in the basis by taking product over all diseases
	DT.no.affy[,q_i:=1-exp(sum(lp0)),by=id]
	## our empirical prior for these diseases is then
	emp<-mean(DT.no.affy[,list(mean_qi=mean(q_i)),by=id]$mean_qi)
	## prior odds
	po<-emp/(1-emp)
	## note that emp us for h1 that beta != 0 therefore we need to take reciprocal as equation assumes pi_0 - see notes
	po<-1/po
	DT.no.affy[,uABF:=po*(q_i/(1-q_i))]
	## add back in affy once we have q values
	DT<-rbind(DT.no.affy,DT[DT$disease == 'aff.t1d',],fill=TRUE)
	DT[,uABF:=max(uABF,na.rm=TRUE),by=id]

	## unfortunately this produces some very large BF somw of which are infinite. One way around this is to set really large and infinite BF
	## to something tractable.
	## I looked at numeric calculations and even with sticking to the log scale calculations will be infinite. For time being set uABF to max(uABF). This will
	## have the effect of averaging over SNPs in a region (even when one is a clear winner) but we may want this to offset any errors in our OR estimation
	## take the strategy that top 0.0001 are actually equal and assign this
	BIG<-quantile(DT[is.finite(DT$uABF),]$uABF,prob=0.9999)
	## BIG is still rather large but at least tractable
	DT[is.infinite(uABF) | uABF > BIG, uABF:=BIG]
	DT[,pwi:=computePWI(uABF,emp),by=c('disease','ld.block')]
	DT[,gh_ss_pw:=gh_ss * pwi]
	DT[,gh_maf_pw:=gh_maf * pwi]

	#DT<-rbind(DT,createControl(DT))

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
cidx<-which(names(bases) %in% c('gh_ss_pp','gh_maf_pp','gh_ss_pw','gh_maf_pw'))
unscaled.comparisons<-lapply(cidx,function(i){
	## need to replace aff.t1d with unscaled version
	b<-bases[[i]]
	bname<-sub('\\_p[pw]$','',names(bases)[i])
	print(names(bases[i]))
	#t1d.gh<-DT[DT$disease=='aff.t1d',c(bname,'id'),with=FALSE]
	idx<-which(rownames(b)=='aff.t1d')
	## get lors and make sure in the correct order
	b[idx,]<-bases[[bname]][idx,]
	#b[idx,]<-t1d.gh[match(t1d.gh$id,colnames(b)),][[bname]]
	comp<-t1d.compare(b)
        with(comp,compareVarWEuc(pca,proj))
})
names(unscaled.comparisons)<-paste('unscaled',names(bases)[cidx],sep='_')

all.comparisons<-c(scaled.comparisons,unscaled.comparisons)
title<-list(expression(hat(beta)),expression(Z),
expression(hat(gamma[ss])),expression(hat(gamma[ss_ppi])),expression(hat(gamma[ss_ppw])),
expression(hat(gamma[maf])),expression(hat(gamma[maf_ppi])),expression(hat(gamma[maf_ppw])),
expression("Non PPi Input Projection"~hat(gamma[ss_ppi])),expression("Non PPi Input Projection"~hat(gamma[ss_ppw])),
expression("Non PPw Input Projection"~hat(gamma[maf_ppi])),expression("Non PPw Input Projection"~hat(gamma[maf_ppw])))
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
#pdf("/home/ob219/git/as_basis/pdf/analysis1.pdf")
multiplot(plotlist=plots,cols=4)
#dev.off()
