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
pi_1<-1e-4

SIM.FILE<-'/home/ob219/scratch/as_basis/figure_data/analysis1_ill_vs_aff_ci.RDS'
SIM.FILE.NULL<-'/home/ob219/scratch/as_basis/figure_data/analysis1_ill_vs_aff_ci_under_null.RDS'
## this takes a long time and should be hived off
if(file.exists(SIM.FILE)){
  empirical_confidence_intervals<-readRDS(SIM.FILE)
}else{
  source("./simulate_aff_ill.R")
}

if(file.exists(SIM.FILE.NULL)){
  empirical_confidence_intervals_null<-readRDS(SIM.FILE.NULL)
}else{
  source("./simulate_aff_ill_null.R")
}

## next load simulations under the null


(load("/scratch/ob219/as_basis/figure_data/analysis1.RData"))
scaled.comparisons<-lapply(bases,function(b){
        comp<-t1d.compare(b)
        with(comp,compareVarWEuc(pca,proj))
})
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
title<-list(lor=expression(hat(beta)),Z=expression(Z),
gh_ss=expression(hat(gamma[ss])),gh_ss_pp=expression(hat(gamma[ss_ppi])),gh_ss_pw=expression(hat(gamma[ss_ppw])),
gh_maf=expression(hat(gamma[maf])),gh_maf_pp=expression(hat(gamma[maf_ppi])),gh_maf_pw=expression(hat(gamma[maf_ppw])),
unscaled_gh_ss_pp=expression("Unscaled"~hat(gamma[ss_ppi])),unscaled_gh_ss_pw=expression("Unscaled"~hat(gamma[ss_ppw])),
unscaled_gh_maf_pp=expression("Unscaled"~hat(gamma[maf_ppi])),unscaled_gh_maf_pw=expression("Unscaled"~hat(gamma[maf_ppw])))
plots<-lapply(names(all.comparisons),function(n){
				mci<-empirical_confidence_intervals[[n]]
        dat<-all.comparisons[[n]]
        dat<-dat[dat!=0]
        dat<-data.frame(disease=names(dat),distance=dat)
        ## hilight ill.t1d and control sets
        dat$hilight<-logical(length=nrow(dat))
        dat[dat$disease %in% c('ill.t1d','control'),]$hilight<-TRUE
        dat$disease<-factor(dat$disease,levels=dat[order(dat$distance),]$disease)
        ggplot(dat,aes(x=disease,y=distance,fill=hilight,color=hilight)) +
				geom_bar(stat="identity") + theme_bw() + xlab("Trait") + ylab("Distance") +
				ggtitle(title[[n]]) + theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
				scale_fill_manual(name="Bar",values=c('white','firebrick'),guide=FALSE) +
				scale_color_manual(name="Bar",values=c('black','black'),guide=FALSE) +
				geom_rect(xmin=0,xmax=length(levels(dat$disease)) + 0.5,ymin=mci[1],ymax=mci[3],alpha=0.1,color=NA,fill='grey90') +
				geom_hline(yintercept=mci[2],color='firebrick')
})

#pdf("/home/ob219/git/as_basis/pdf/analysis1_w_sim.pdf",width=17,height=12)
#pdf("/home/ob219/git/as_basis/pdf/analysis1_w_sim.pdf",paper="a4")
multiplot(plotlist=plots,cols=3)
# wysiwyg printing adjust on screen and then run
dev.print(pdf,"/home/ob219/git/as_basis/pdf/analysis1_aff_vs_ill_w_sim.pdf")
#dev.off()
