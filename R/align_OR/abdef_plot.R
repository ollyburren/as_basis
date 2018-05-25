## decrufted version of egpa plot where basis
## is computed using ws_emp

library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)

SHRINKAGE_METHOD<-'ws_emp'
## snpsStat1 1kg that have been filtered tio include just basis SNPs
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_feb/'

## prep the data

DT.abdef <- fread("/home/ob219/rds/hpc-work/pid/GWAS_META/AbDeficiency2.chr1-22.assoc.logistic.plink")
#DT.li <- fread("/home/ob219/tmp/CVID_Ichip/CVID_QCed_assoc.logistic")
support.dir<-'/home/ob219/rds/hpc-work/as_basis/support_tab'
ref_af_file<-file.path(support.dir,'as_basis_snp_support_feb_2018_w_ms.tab')
m.DT <- fread(ref_af_file)

DT.abdef[,pid:=paste(CHR,BP,sep=':')]
DT.abdef.filt <- DT.abdef[pid %in% m.DT$pid,]
setkey(DT.abdef.filt,pid)
setkey(m.DT,pid)
## check alleles
m<-DT.abdef.filt[m.DT][!is.na(STAT),]
m[ref_a2==A1,flip:=FALSE]
m[ref_a1==A1,flip:=TRUE]
m[flip==TRUE,OR:=1/OR]

## what we want to project

abdef.dat<-m[,.(or=OR,trait='abdef',pid)]
ref.DT <- m.DT[pid %in% abdef.dat$pid,]
write.table(ref.DT,file="/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_abdef.tab",row.names=FALSE,quote=FALSE,sep="\t")

## process the dataset so that we can project


## import GWAS data for basis
## support files

# reference allele frequencies
ref_af_file<-"/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_abdef.tab"
#ld_file<-file.path(support.dir,'all.1cM.tab')
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
## project on biobank to see if we can recreate Chris' figure.
mat.emp <- create_ds_matrix(abdef.dat,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=mat.emp)
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

analytical.vars <- compute_proj_var(gw.DT,w.DT,shrink.DT,DEFAULT.SNPSTATS.DIR,SHRINKAGE_METHOD,quiet=TRUE)
## we can convert between by normalising by our factor total.num/(ncases.num * (total.num-ncases.num))
N <- gw.DT$n[1]
N1 <- gw.DT$n1[1]
factor <- N/(N1 * (N-N1))
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars/factor)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)

## integrate with projections

mani <- fread(m_file)[,.(trait,cases,controls,basis_trait)]
mani <- rbind(mani,data.table(trait="abdef",cases=775,controls=9958-775,basis_trait=FALSE))
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
#merged[,glabel:=gsub("egpa\\_","",trait)]
setnames(merged,'variable','pc')

## for biplot we need to dcast

biplot.DT <- dcast(merged,glabel~pc)
## need to add in CI's also
ci.DT <- melt(merged,id.vars=c('glabel','pc'),measure.vars=c('ci_lo','ci_hi')) %>% dcast(.,glabel~pc+variable)
biplot.DT <- cbind(biplot.DT,ci.DT)



library(cowplot)
library(ggrepel)
PC1.var<-signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3)
PC2.var<-signif(summary(pc.emp)[['importance']][2,]["PC2"]*100,digits=3)
pp<-ggplot(biplot.DT,aes(x=PC1,y=PC2,label=glabel,ymin=PC2_ci_lo,ymax=PC2_ci_hi,xmin=PC1_ci_lo,xmax=PC1_ci_hi,col=glabel=="abdef")) +
geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5) +
geom_point(size=3) + geom_text_repel()   + scale_color_manual(guide=FALSE,values=c('black','firebrick')) +
scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + coord_cartesian(xlim=c(-0.15,0.16)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +
background_grid(major = "xy", minor = "none")

save_plot("~/tmp/abdef_bi.pdf",pp)

## non biplot

trait.hilight <- c('anca_Neg','egpa','mpo','mpo_Pos')
#trait.hilight <- c('jia_ERA','mpo_Pos')

## label max an min by pc
mins <- merged[,.SD[which.min(value),],by='pc']
maxs <- merged[,.SD[which.max(value),],by='pc']
merged[value %in% c(mins$value,maxs$value),label:=trait]
merged[is.na(label),label:='']
merged[trait %in% trait.hilight,cat:=trait]

cDT<-merged[trait=='control',.(pc,control.v=value)]
setkey(merged,'pc')
setkey(cDT,'pc')
merged<-cDT[merged]
#DT[,value:=value-control.v]
merged[,pc.plot:=factor(pc,levels=paste0('PC',1:11))]

ppa <- ggplot(merged,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + scale_color_discrete ("Disease Subtype",labels=c('ANCA -','EGPA','MPO','MPO +','Other')) +
ggtitle("Variance calculated under null")

## we stop here because if we empirically estimate the variance then it is huge, this is a consequence of poor genomic control in the original study

stop()

if(FALSE){
  #plot without substracting control !
ggplot(merged,aes(x=pc,y=value,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ci_lo,ymax=ci_hi)) +
geom_line(position = position_dodge(width = 0.9)) + geom_text_repel() + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none")
}
#save_plot("~/tmp/egpa_plot_with_ci.pdf",ppg,base_height=5)

## what happens if we analytically compute the variance using actual standard errors
if(FALSE){
  ## takes a while so do once and save
  by.trait <- split(egpa.DT,egpa.DT$trait)
  DT <- by.trait[[1]]
  ## compute standard error from OR and p.value
  aly.vars <- lapply(by.trait,function(DT){
    message(sprintf("Processing %s",unique(DT$trait)))
    DT[,Z:=qnorm(p.value/2,lower.tail=FALSE)]
    DT[,emp_se:=log(or)/Z]
    ## for those that are uncomputable estimate
    DT[is.na(emp_se) | is.infinite(emp_se),emp_se:=se_null(n,n1,maf)]
    DT[,c('chr','position'):=tstrsplit(pid,':')]
    tmp<-compute_proj_var(DT,w.DT,shrink.DT,DEFAULT.SNPSTATS.DIR,SHRINKAGE_METHOD,quiet=TRUE)
    data.table(trait=unique(DT$trait),pc=names(tmp),vars=tmp)
  })

  saveRDS(aly.vars,file="~/tmp/alyegpa.vars.RDS")
}
aly.vars <- readRDS("~/tmp/alyegpa.vars.RDS") %>% rbindlist
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
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + scale_color_discrete ("Disease Subtype",labels=c('ANCA -','EGPA','MPO','MPO +','Other')) +
ggtitle("Variance calculated under alternative")

## do locally so quicker
plot_grid(ppa, ppb,labels = c("A", "B"),nrow=2)
## consolidated without the titles for thesis
ppa <- ggplot(merged,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + #scale_color_discrete ("Disease Subtype",labels=c('ANCA -','EGPA','MPO','MPO +','Other')) +
theme(legend.position="bottom") + coord_cartesian(ylim=c(-.2,.2))

ppb <- ggplot(merged2,aes(x=pc.plot,y=value-control.v,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=ifelse(ci_lo==0,0,ci_lo-control.v),ymax=ifelse(ci_hi==0,0,ci_hi-control.v))) +
geom_line(position = position_dodge(width = 0.9)) + geom_text(show.legend = FALSE) + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none") + xlab("PC") + ylab("PC Loading - Control") + guides(alpha=FALSE) + #scale_color_discrete ("Disease Subtype",labels=c('ANCA -','EGPA','MPO','MPO +','Other')) +
theme(legend.position="bottom") + coord_cartesian(ylim=c(-.2,.2))

p <- plot_grid(ppa, ppb,labels = c("A", "B"),nrow=2)

save_plot("~/tmp/egpa_all_pc_with_ci.pdf",p,base_width=8,base_height=8)



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
  M <- M[trait %in% c(mani[basis_trait!=1,]$trait),]
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
  pdf("~/tmp/pid_tree_dist_z.pdf",height=7,width=7)
  plot(hc, hang=-1,xlab="", sub="",cex=1.1,main="Z score")
  dev.off()
}else{
  pdf("~/tmp/pid_tree.pdf",height=7,width=7)
  plot(hc, hang=-1,xlab="", sub="",cex=1.1,main="Raw")
  dev.off()
}
