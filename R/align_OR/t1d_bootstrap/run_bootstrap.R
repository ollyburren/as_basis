## run_bootstrap projection
library(optparse)

TEST <- FALSE

## we expect that the bs_matrix file is in the same dir as relevant support files
## these need to be called samples.RDS and snps.RDS

## THESE ARE HARDCODED FOR NOW
BASIS.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/basis_ws.RDS"
SHRINK.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/shrink_ws.RDS"
SUPPORT.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/snp_support_bootstrap_USE.tab"

option_list = list(
        make_option(c("-m", "--bs_matrix"), type="character", default=NULL,
              help="RDS file containing sample bootstraps", metavar="character")
      )


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if(!file.exists(args$bs_matrix))
    stop("Can't find source data files")
}else{
  #for testing
  message("IN TESTING MODE ======>!")
  args <- list(
    bs_matrix='/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims/wtccc_1000_1000_500/bs_34.RDS'
  )
}


library(snpStats)
library(cupcake)

bs.idx <- readRDS(file.path(args$bs_matrix))
IN.DIR<-dirname(args$bs_matrix)
# read in bs index
samples <- readRDS(file.path(IN.DIR,'samples.RDS')) %>% as.data.frame
rownames(samples) <- samples$sampleid
snps <- readRDS(file.path(IN.DIR,'snps.RDS'))
#get pca for projections
basis <- readRDS(BASIS.FILE)
shrink <- readRDS(SHRINK.FILE)
support <- fread(SUPPORT.FILE)
setkey(support,pid)
#conduct logistic fit
res.dt <- lapply(1:ncol(bs.idx),function(bs.col){
#res.dt <- lapply(1:2,function(bs.col){
  message(sprintf("Running GWAS %d",bs.col))
  idx <- bs.idx[,bs.col]
  bs.samples <- samples[idx,]
  bs.snps <- snps[idx,]
  n1 <- sum(bs.samples$cc==1)
  n <- nrow(bs.samples)
  rownames(bs.snps) <- rownames(bs.samples)
  res <- snp.rhs.estimates(cc ~ b58region,family='binomial',data=bs.samples,snp.data=bs.snps,uncertain = TRUE)
  or <- sapply(res,'[[','beta')
  p <- sapply(res,function(x) 2*pnorm(abs(x$beta/sqrt(x$Var.beta)),lower.tail=FALSE))
  z<-sapply(res,function(x) x$beta/sqrt(x$Var.beta))
  tmp.dt <- data.table(pid=names(res),or=exp(or),p.value=p,trait=paste("bs",bs.col,sep='_'),n=n,n1=n1)
  setkey(tmp.dt,pid)
  ## next we need to add in annotation data from support file
  support[tmp.dt]
}) %>% rbindlist
setkey(res.dt,pid)
res.mat.emp<-create_ds_matrix(res.dt,shrink,'ws_emp')
pred.emp <- predict(basis,newdata=res.mat.emp)
ofile <- file.path(IN.DIR,gsub("^bs","proj",basename(args$bs_matrix)))
saveRDS(pred.emp,file=ofile)


## create a list of commands to execute !



if(FALSE){
library(cowplot)
library(ggrepel)
library(magrittr)
library(cupcake)
BASIS.FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/basis.RDS"
basis <- readRDS(BASIS.FILE)
id <- '/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims/wtccc_1900_1900_500'
pred.emp<- lapply(list.files(path=id,pattern="proj_*",full.names=TRUE),readRDS) %>% do.call('rbind',.)
rownames(pred.emp) <- paste('t1d',1:nrow(pred.emp),sep='_')
DT<-data.table(trait=rownames(basis$x),PC1=basis$x[,"PC1"],PC2=basis$x[,"PC2"],basis.trait=T)
library(ggplot2)
DT.res <- data.table(trait="",PC1=pred.emp[,"PC1"],PC2=pred.emp[,"PC2"],basis.trait=F)
ggplot(rbind(DT,DT.res),aes(x=PC1,y=PC2,label=trait,color=basis.trait)) + geom_point() +
geom_text_repel() + scale_color_manual(guide=FALSE,values=c('firebrick','dodgerblue'))

samps<-list.files(path='/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims_ws',pattern="wtccc*",full.names=TRUE)

all.res.alt <- lapply(samps,function(path){
  message(path)
  vars <- lapply(list.files(path=path,pattern="^proj*",full.names=TRUE),readRDS) %>% do.call('rbind',.) %>% apply(.,2,var)
  cases <- gsub("wtccc\\_([^_]+)\\_(.*)","\\1",basename(path))
  controls <- gsub("wtccc\\_([^_]+)\\_([^_]+)\\_.*","\\2",basename(path))
  data.table(cases=cases,controls=controls,vars,pc=names(vars))
}) %>% rbindlist



OUTDIR <- '/home/ob219/rds/hpc-work/as_basis/callibration/filter_bootstrap_w_ms_ws_fixed'
af <- list.files(path=OUTDIR,pattern="*.RDS",full.names=TRUE)
byexp<-split(af,gsub("\\_[a-z]+.RDS$","",basename(af)))


## build callibration curve for one sample size to start to see if there is a trend

all.res <- lapply(names(byexp),function(x){
  all.res <- lapply(byexp[[x]],readRDS) %>% do.call('rbind',.)
  var <- apply(all.res,2,var)
  tmp<-strsplit(x,"\\_") %>% unlist
  data.table(total=tmp[1],ncases=tmp[2],variance=var,pc=names(var))
}) %>% rbindlist

all.res <- all.res[order(as.numeric(total)),]

## it does not make sense to have 1900 cases and 100 controls so delete !
which.idx <- which(all.res$total==2000 & all.res$ncases==1900)
all.res <- all.res[-which.idx,]




## if our estimate of the variance is drawn
CQF.lt <- qchisq(0.025, 500-1, lower.tail=FALSE)
CQF.ut <- qchisq(0.025, 500-1, lower.tail=TRUE)

## an approximate way to do this is variance * sqrt(2/(n-1)) ## I should ask Chris how this works



all.res[,c('ci.lower','ci.upper'):=list((variance * (500-1)/CQF.lt),((variance * (500-1))/CQF.ut)) ]

all.res[,ci:=variance * sqrt(2/(500-1))]
all.res[,total:=factor(total,levels=sort(unique(as.numeric(total))))]

## plot for PC1

# pc1 <- all.res[pc=='PC1',list(m.var=mean(variance)),by=c('total')]
# pc1[,c('ci.lower','ci.upper'):=list((m.var * (500-1)/CQF.lt),((m.var * (500-1))/CQF.ut)) ]
# pc1[,total:=as.numeric(as.character(total))]

all.res.alt[,total:=as.numeric(cases)+as.numeric(controls)]
all.res.alt<-all.res.alt[order(total),]

pc1.alt <- all.res.alt[pc=='PC1',]
pc1.alt[,total:=as.numeric(cases)+as.numeric(controls)]

library(ggplot2)
library(scales)
library(cowplot)

#all.res[,c('ci.lower','ci.upper'):=list((variance * (500-1)/CQF.lt),((variance * (500-1))/CQF.ut)) ]

#all.res[,ci:=variance * sqrt(2/(500-1))]
#all.res[,total:=as.numeric(total)]
#all.res[,ncases:=as.numeric(ncases)]
#all.res[,total:=factor(total,levels=sort(unique(as.numeric(total))))]
pd <- position_dodge(width=0.1)
ppf <- ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + xlab("Cases") + scale_color_discrete("Sample Size") + ylab("Variance (PC1 Loadings)") +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + theme(legend.position =c(0.65,0.75)) +
background_grid(major = "xy", minor = "none")

#ggplot(pc1,aes(x=log(total),y=log(m.var),ymin=log(m.var-ci.lower), ymax=log(m.var+ci.upper))) + geom_errorbar(width=.1) +
#geom_point(col='dodgerblue') + geom_smooth(method = "glm", formula=y ~ x) + geom_point(data=pc1.alt,aes(x=log(total),y=log(vars)),col='firebrick',inherit.aes=FALSE)

#ggplot(pc1,aes(x=total,y=m.var,ymin=m.var-ci.lower, ymax=m.var+ci.upper)) + geom_errorbar(width=.1) +
#geom_point(col='dodgerblue') + geom_smooth(method = "glm", formula=y ~ x) + geom_point(data=pc1.alt,aes(x=total,y=vars),col='firebrick',inherit.aes=FALSE) +
#scale_x_continuous(trans="log",breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x))) +
#scale_y_continuous(trans="log",breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x))) + xlab("Sample Size") +
#ylab("Variance (PC1 Loading)")



ppg <- ggplot(pc1[total!=350000,],aes(x=total,y=m.var,ymin=m.var-ci.lower, ymax=m.var+ci.upper)) + geom_errorbar(width=.1) +
geom_point(col='dodgerblue') + geom_smooth(method = "glm", formula=y ~ x) + geom_point(data=pc1.alt,aes(x=total,y=vars),col='firebrick',inherit.aes=FALSE) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))  +
xlab("Sample Size") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none")

plots <- plot_grid(ppf, ppg, labels = c("A", "B"))

save_plot("~/tmp/var_simulation_null.pdf",plots,base_width = 10)

ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + geom_point(data=all.res.alt[pc=='PC1',],aes(x=total,y=vars),col='firebrick',inherit.aes=FALSE) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
xlab("Sample Size") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none")


## what about variance scaling how well does it fit the bootstrap ?
all.res.alt[,total:=factor(total,levels=levels(all.res$total))]
all.res.alt[,total.num:=as.numeric(as.character(total))]
all.res.alt[,ncases.num:=as.numeric(cases)]
all.res.alt[,factor:=total.num/(ncases.num * (total.num-ncases.num))]

all.res[,total.num:=as.numeric(as.character(total))]
all.res[,ncases.num:=as.numeric(ncases)]
all.res[,factor:=total.num/(ncases.num * (total.num-ncases.num))]


#all.vars[,prop:=ncases/as.numeric(as.character(total))]
#all.vars[,factor:=prop * (1-prop)]

## derive the variance for all others using 350000
r.DT <- subset(all.res, total==350000 & ncases==50000 & pc=='PC1')
rf.factor <- r.DT$factor
rf.variance <- r.DT$variance
all.res.alt[,dvar:=(factor/rf.factor) * rf.variance]

## load in the actual variances that we compute in

ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) +
geom_point(data=all.res.alt[pc=='PC1',],aes(x=as.numeric(cases),y=vars,color=total),inherit.aes=FALSE,size=3) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
xlab("# Cases") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none") +
geom_point(data=all.res.alt[pc=='PC1',],aes(x=as.numeric(cases),y=dvar,color=total),pch=17,size=5,alpha=0.3,inherit.aes=FALSE)


## next load in the empirical variances that we calculate with run_actual_gwas_withSE_variance_calculation

var.DT <- readRDS("~/tmp/act_var_bootstrap.RDS")

## we can make nicer plots if we consolidate all.res.alt (both types) and var.data

setnames(var.DT,c('cases','total','pc','variance'))
var.DT[,ptype:='GWAS_SE_calculated']

## next we format bootstrap_empirical

bs.emp <- all.res.alt[,.(cases,total,pc,vars)]
setnames(bs.emp,'vars','variance')
bs.emp[,ptype:='bootstrap_empirical']

## next we format bootstrap estimated from null

bs.emp.null <- all.res.alt[,.(cases,total,pc,dvar)]
setnames(bs.emp.null,'dvar','variance')
bs.emp.null[,ptype:='null_se_calculated']

vcals.DT <- rbindlist(list(var.DT,bs.emp,bs.emp.null))



var.DT[,total:=factor(total,levels=levels(all.res$total))]
ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) +
geom_point(data=all.res.alt[pc=='PC1',],aes(x=as.numeric(cases),y=vars,color=total),pch=15,inherit.aes=FALSE,size=3) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
xlab("# Cases") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none") +
geom_point(data=var.DT[variable=='PC1',],aes(x=as.numeric(cases),y=value,color=total),pch=17,size=5,alpha=0.3,inherit.aes=FALSE,position = position_dodge(width = 0.9)) +
geom_point(data=all.res.alt[pc=='PC1',],aes(x=as.numeric(cases),y=dvar,color=total),pch=8,size=5,inherit.aes=FALSE,position = position_dodge(width = 0.9))

## remove confs that are not in the bootstrap to make plotting clearer



ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
xlab("# Cases") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none") +
geom_point(data=vcals.DT[pc=='PC1' & total!=350000],aes(x=as.numeric(cases),y=variance,color=total,pch=ptype),size=3,position=pd,inherit.aes=FALSE)



## what about other PC's

pcs <- all.res[,list(m.var=mean(variance)),by=c('total','pc')]
pcs[,c('ci.lower','ci.upper'):=list((m.var * (500-1)/CQF.lt),((m.var * (500-1))/CQF.ut)) ]
pcs[,total:=as.numeric(as.character(total))]

pcs.alt <- all.res.alt
pcs.alt[,total:=as.numeric(cases)+as.numeric(controls)]

baral <- 0.1

pcs[,pc:=factor(pc,levels=paste0('PC',1:11))]
fpg <- ggplot(pcs[total!=350000 & !pc %in% c('PC10','PC11'),],aes(x=total,y=m.var,ymin=m.var-ci.lower, ymax=m.var+ci.upper)) + geom_errorbar(width=.1,alpha=baral,col='dodgerblue') +
geom_point(col="dodgerblue") + geom_line(stat="smooth",method = "glm", formula=y ~ x,se = FALSE,alpha=0.3,col='dodgerblue') + geom_point(data=pcs.alt[pc!='PC10',],aes(x=total,y=vars),col='firebrick',inherit.aes=FALSE) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))  +
xlab("Sample Size") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none") + facet_wrap(~pc)
save_plot("~/tmp/all_pcs_calli.pdf",fpg,base_height=10)

## obtain the callibration curve by fitting to all data apart from the last pc

dat <- all.res[pc!='PC11',]

dat[,c('l.ss','l.var'):=list(log(as.numeric(as.character(total))),log(variance))]

calli <- lm(l.var ~ l.ss,data=dat)

saveRDS(calli,file="/home/ob219/rds/hpc-work/as_basis/callibration/callibration_model.RDS")


}
