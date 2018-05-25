library(cowplot)
library(magrittr)
library(cupcake)
library(ggplot2)
library(scales)



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



pd <- position_dodge(width=0.1)
## what about variance scaling how well does it fit the bootstrap ?
all.res.alt[,total:=factor(total,levels=levels(all.res$total))]
all.res.alt[,total.num:=as.numeric(as.character(total))]
all.res.alt[,ncases.num:=as.numeric(cases)]
all.res.alt[,factor:=total.num/(ncases.num * (total.num-ncases.num))]

all.res[,total.num:=as.numeric(as.character(total))]
all.res[,ncases.num:=as.numeric(ncases)]
all.res[,factor:=total.num/(ncases.num * (total.num-ncases.num))]

r.DT <- subset(all.res, total==350000 & ncases==50000 & pc=='PC1')
rf.factor <- r.DT$factor
rf.variance <- r.DT$variance
all.res.alt[,dvar:=(factor/rf.factor) * rf.variance]


## next load in the empirical variances that we calculate with run_actual_gwas_withSE_variance_calculation

var.DT <- readRDS("~/tmp/act_var_bootstrap.RDS")

## we can make nicer plots if we consolidate all.res.alt (both types) and var.data

setnames(var.DT,c('cases','total','pc','variance'))
var.DT[,ptype:='GWAS_se_calculated']

## next we format bootstrap_empirical

bs.emp <- all.res.alt[,.(cases,total,pc,vars)]
setnames(bs.emp,'vars','variance')
bs.emp[,ptype:='bootstrap_empirical']

## next we format bootstrap estimated from null

#bs.emp.null <- all.res.alt[,.(cases,total,pc,dvar)]
#setnames(bs.emp.null,'dvar','variance')
#bs.emp.null[,ptype:='null_se_calculated']

all.res[,dvar:=(factor/rf.factor) * rf.variance]

bs.calc.null <- all.res[,.(ncases,total,pc,dvar)]
setnames(bs.calc.null,'dvar','variance')
setnames(bs.calc.null,'ncases','cases')
bs.calc.null[,ptype:='null_se_calculated']



#vcals.DT <- rbindlist(list(var.DT,bs.emp,bs.emp.null))
vcals.DT <- rbindlist(list(var.DT,bs.emp))



## remove confs that are not in the bootstrap to make plotting clearer

#plot locally so quicker !
save(list=c('all.res','vcals.DT','bs.calc.null'),file="~/tmp/simulation_curves.RData")
load("~/tmp/simulation_curves.RData")
pd <- position_dodge(width=0.1)
pp<-ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
#geom_point(position=pd) + geom_line(position=pd,alpha=0.2) + geom_errorbar(width=.1,position=pd,,alpha=0.2) +
geom_point(position=pd) +  geom_errorbar(width=.1,position=pd,,alpha=0.2) +
scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
xlab("# Cases") + ylab("Variance (PC1 Loading)") + background_grid(major = "xy", minor = "none") +
geom_point(data=vcals.DT[pc=='PC1' & total!=350000],aes(x=as.numeric(cases),y=variance,color=total,pch=ptype),size=3,position=position_dodge(width=0.2),inherit.aes=FALSE) +
geom_point(data=bs.calc.null[pc=='PC1' & total!=350000],aes(x=as.numeric(cases),y=variance,color=total),pch=3,size=3,position=pd,inherit.aes=FALSE) +
theme(legend.position="bottom") + scale_color_discrete ("Study Size") + scale_shape_manual(values=c(15,17)) + guides(pch=FALSE)
save_plot(pp,file="~/tmp/var_sims.pdf",base_width=7)
