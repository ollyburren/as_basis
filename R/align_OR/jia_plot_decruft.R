## decrufted version of jia plot where basis
## is computed using ws_emp

library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)

SHRINKAGE_METHOD<-'ws_emp'
## snpsStat1 1kg that have been filtered tio include just basis SNPs
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_feb/'


## import GWAS data for basis
## support files
support.dir<-'/home/ob219/rds/hpc-work/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support_feb_2018_w_ms.tab')
#ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
## project on biobank to see if we can recreate Chris' figure.
jia_traits<-fread(m_file)[grep('jia|pso|bb_RA|bb_AS',trait),]$trait
## add in pso traits
jia.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,TRUE,jia_traits)
jia.mat.emp<-create_ds_matrix(jia.DT,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=jia.mat.emp)
emp<-rbind(pc.emp$x,pred.emp)
emp <- cbind(as.data.table(emp),trait=rownames(emp))


## compute and add confidence intervals
## use null to start with as quicker

## get an exemplar study

gw.DT <- basis.DT[trait==head(sample(unique(trait)),n=1),]
## set the standard error to the null
gw.DT[,emp_se:=se_null(n,n1,maf)]
## data table of PC snp loadings
gw.DT[,c('chr','position'):=tstrsplit(pid,':')]
w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)

analytical.vars <- compute_proj_var(gw.DT,w.DT,shrink.DT,DEFAULT.SNPSTATS.DIR,SHRINKAGE_METHOD,quiet=FALSE)
## we can convert between by normalising by our factor total.num/(ncases.num * (total.num-ncases.num))
N <- gw.DT$n[1]
N1 <- gw.DT$n1[1]
factor <- N/(N1 * (N-N1))
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars/factor)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)

## integrate with projections

mani <- fread(m_file)[,.(trait,cases,controls,basis_trait)]
emp <- melt(emp,id.vars=c('trait'))
setkey(emp,trait)
setkey(mani,trait)
merged <- mani[emp]
## avoid integer overflow
merged[,c('cases','controls'):=list(as.numeric(cases),as.numeric(controls))]
merged[,var.factor:=(cases+controls)/(cases*controls)]
setkey(merged,'variable')
merged <- merged[adj.vars]
merged[,variance_null:=mfactor * var.factor]
## next compute the variance
merged[,ci:=sqrt(variance_null)*1.96]
merged[,c('ci_lo','ci_hi'):=list(value - ci,value + ci)]
merged[basis_trait==1,c('ci_lo','ci_hi'):=list(NA,NA)]

merged[,glabel:=trait]
merged[,glabel:=gsub("jia\\_","",trait)]
setnames(merged,'variable','pc')

## for biplot we need to dcast

biplot.DT <- dcast(merged,glabel~pc)
## need to add in CI's also
ci.DT <- melt(merged,id.vars=c('glabel','pc'),measure.vars=c('ci_lo','ci_hi')) %>% dcast(.,glabel~pc+variable)
setnames(ci.DT,"glabel","trait")
biplot.DT <- cbind(biplot.DT,ci.DT[])

## get a list of basis traits and use this to set a label

biplot.DT[trait %in% mani[basis_trait==1,]$trait,category:=glabel]


library(cowplot)
library(ggrepel)
PC1.var<-signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3)
PC2.var<-signif(summary(pc.emp)[['importance']][2,]["PC2"]*100,digits=3)
pp <- ggplot(biplot.DT,aes(x=PC1,y=PC2,label=glabel,ymin=PC2_ci_lo,ymax=PC2_ci_hi,xmin=PC1_ci_lo,xmax=PC1_ci_hi,col=category,alpha=is.na(category))) +
geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5) +
#geom_point(size=3) + geom_text_repel()   + scale_color_manual(guide=FALSE,values=c('black','grey')) +
geom_point(size=3) + geom_text_repel() +
scale_alpha_discrete(guide=FALSE,range=c(0.5,1)) + coord_cartesian(xlim=c(-0.15,0.16)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +
background_grid(major = "xy", minor = "none") + guides(col=FALSE)
save_plot("~/tmp/jia_biplot.pdf",pp,base_aspect_ratio=1,base_width=7,base_height=7)

## next do the hclust
zDist <- FALSE
controlv <- merged[trait=='control',.(pc,covalue=value)]
td <- merged[,.(trait,pc,value,ci)]
setkey(td,'pc')
setkey(controlv,'pc')
M <- controlv[td]
## add variance explained
vexp<-summary(pc.emp)$importance[2,]
vexp <- data.table(pc=names(vexp),vexp=vexp)
setkey(vexp,'pc')
M <- vexp[M]
if(zDist){
  M[trait=='control',ci:=1]
  #M[,dist:=(value-covalue)]
  M[,dist:=(value-covalue)/(ci/1.96)]
  M <- M[trait %in% c(mani[basis_trait!=1,]$trait,'control'),]
}else{
  M[,dist:=(value-covalue)]
}
dist<-melt(M,id.vars=c('trait','pc'),measure.vars='dist') %>% dcast(.,trait ~ pc)
foo <- as.matrix(dist[,-1])
rownames(foo)<-dist$trait
hc <- hclust(dist(foo))
cols <- rep("black",length(hc$labels))
cols[grep("^bb_",hc$labels)]<-'olivedrab4'

par(mar=c(0, 4, 4, 2)) # c(bottom, left, top, right)
if(zDist){
  pdf("~/tmp/jia_tree_dist_z.pdf",height=7,width=7)
  plot(hc, hang=-1,xlab="", sub="",cex=1.1,main="Z score")
  dev.off()
}else{
  pdf("~/tmp/jia_tree.pdf",height=7,width=7)
  plot(hc, hang=-1,xlab="", sub="",cex=1.1,main="Raw")
  dev.off()
}

## non biplot

## label max an min by pc
mins <- merged[,.SD[which.min(value),],by='pc']
maxs <- merged[,.SD[which.max(value),],by='pc']
merged[value %in% c(mins$value,maxs$value),label:=trait]
merged[is.na(label),label:='']
merged[grep('jia',trait),cat:=trait]

cDT<-merged[trait=='control',.(pc,control.v=value)]
setkey(merged,'pc')
setkey(cDT,'pc')
merged<-cDT[merged]
#DT[,value:=value-control.v]
merged[,pc.plot:=factor(pc,levels=paste0('PC',1:11))]

ppa <- ggplot(merged,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + scale_color_discrete ("Disease Subtype",labels=c('ERO','ERA','PO','PsA','RF-','RF+','Sys','Other')) +
ggtitle("Variance calculated under null")

if(FALSE){
  #plot without substracting control !
ggplot(merged,aes(x=pc,y=value,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ci_lo,ymax=ci_hi)) +
geom_line(position = position_dodge(width = 0.9)) + geom_text_repel() + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none")
}
#save_plot("~/tmp/jia_plot_with_ci.pdf",ppg,base_height=5)

## what happens if we analytically compute the variance using actual standard errors
if(FALSE){
  ## takes a while so do once and save
  by.trait <- split(jia.DT,jia.DT$trait)
  DT <- by.trait[[1]]
  ## compute standard error from OR and p.value
  aly.vars <- lapply(by.trait,function(DT){
    message(sprintf("Processing %s",unique(DT$trait)))
    DT[,Z:=qnorm(p.value/2,lower.tail=FALSE)]
    DT[,emp_se:=log(or)/Z]
    ## for those that are uncomputable estimate
    DT[is.na(emp_se),emp_se:=se_null(n,n1,maf)]
    DT[,c('chr','position'):=tstrsplit(pid,':')]
    tmp<-compute_proj_var(DT,w.DT,shrink.DT,DEFAULT.SNPSTATS.DIR,SHRINKAGE_METHOD,quiet=TRUE)
    data.table(trait=unique(DT$trait),pc=names(tmp),vars=tmp)
  })

  saveRDS(aly.vars,file="~/tmp/aly.vars.RDS")
}
aly.vars <- readRDS("~/tmp/aly.vars.RDS") %>% rbindlist
setkeyv(aly.vars,c('trait','pc'))
merged2 <- copy(merged)
setkeyv(merged2,c('trait','pc'))

merged2<-aly.vars[merged2]
## recompute the ci etc.
merged2[,ci:=sqrt(vars)*1.96]
merged2[,c('ci_lo','ci_hi'):=list(value - ci,value + ci)]
merged2[basis_trait==1,c('ci_lo','ci_hi'):=list(NA,NA)]

save(list=c('merged','merged2'),file="~/tmp/plott_pc.RData")

ppb <- ggplot(merged2,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + scale_color_discrete ("Disease Subtype",labels=c('ERO','ERA','PO','PsA','RF-','RF+','Sys','Other')) +
ggtitle("Variance calculated under alternative")

## do locally so quicker
plot_grid(ppa, ppb,labels = c("A", "B"),nrow=2)
## consolidated without the titles for thesis
ppa <- ggplot(merged,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + scale_color_discrete ("Disease Subtype",labels=c('ERO','ERA','PO','PsA','RF-','RF+','Sys','Other')) +
theme(legend.position="bottom")

ppb <- ggplot(merged2,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + scale_color_discrete ("Disease Subtype",labels=c('ERO','ERA','PO','PsA','RF-','RF+','Sys','Other')) +
theme(legend.position="bottom")

p <- plot_grid(ppa, ppb,labels = c("A", "B"),nrow=2)

save_plot("~/tmp/jia_all_pc_with_ci.pdf",p,base_width=8,base_height=8)
