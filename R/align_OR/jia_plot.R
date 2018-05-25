library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)

SHRINKAGE_METHOD<-'emp'
#SHRINKAGE_METHOD<-'ws_emp'
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
jia_traits<-fread(m_file)[grep('jia|pso',trait),]$trait
## add in pso traits
jia.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,TRUE,jia_traits)
jia.mat.emp<-create_ds_matrix(jia.DT,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=jia.mat.emp)
emp<-rbind(pc.emp$x,pred.emp)
ml<-list(
  CD = 'bb_CD',
  CEL = 'bb_CEL',
  MS = 'bb_MS',
  RA = 'bb_RA',
  SLE = 'bb_SLE',
  T1D = 'bb_T1D',
  UC = 'bb_UC',
  PBC = 'PBC',
  PSC = 'PSC',
  asthma = 'asthma'
)
g <- function(M){
    M <- cbind(as.data.table(M),trait=rownames(M))
    M$compare<-"none"
    for(i in seq_along(ml)) {
        M[trait %in% c(names(ml)[i], ml[i]), compare:=names(ml)[i]]
    }
    M[trait=="control",compare:="control"]
    M
}
emp<-g(emp)[,is.biobank:=FALSE]
emp[grep("^bb",trait),is.biobank:=TRUE]
emp[,is.imb.ai:=compare %in% names(ml)]
emp[,glabel:=trait]
emp[is.imb.ai==FALSE,glabel:='']
emp[trait=='control',glabel:='CONTROL']
emp[,glabel:=gsub("jia\\_","",trait)]

## add confidence intervals

melt(emp,id.vars=c('trait','total','ncases'),measure.vars=names(emp)[grep('PC',names(emp))])

## get sample sizes from manifest

mani <- fread(m_file)[,.(trait,cases,controls)]
mani <- mani[,c('total','ncases'):=list(cases+controls,cases)]
## load in model for variance callibration
callibration <- readRDS("/home/ob219/rds/hpc-work/as_basis/callibration/callibration_model_fixed.RDS")

mani[,basis.var:=exp(predict(callibration,mani))]
setkey(mani,trait)
setkey(emp,trait)

merged <- mani[emp]
merged[trait=='control',basis.var:=0]
merged[,PC.ci:=sqrt(basis.var)*1.96]


ggplot(emp,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Empirical MAF SE shrinkage')
## what happens if we use estimate instead ?
library(cowplot)
library(ggrepel)
PC1.var<-signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3)
PC2.var<-signif(summary(pc.emp)[['importance']][2,]["PC2"]*100,digits=3)
ppf<-ggplot(emp,aes(x=PC1,y=PC2,color=is.imb.ai,label=glabel)) + geom_point(size=3) + geom_text_repel()   +
scale_color_manual(guide=FALSE,values=c('black','grey')) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + coord_cartesian(xlim=c(-0.15,0.16)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none")
save_plot("~/tmp/jia_plot.pdf",ppf)

ppg<-ggplot(merged,aes(x=PC1,y=PC2,color=is.imb.ai,label=glabel,ymin=PC2-PC.ci,ymax=PC2+PC.ci,xmin=PC1-PC.ci,xmax=PC1+PC.ci)) + geom_point(size=3) + geom_text_repel()   +
scale_color_manual(guide=FALSE,values=c('black','grey')) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + coord_cartesian(xlim=c(-0.15,0.16)) +
xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none") + geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5)
save_plot("~/tmp/jia_plot_with_ci.pdf",ppg,base_height=5)


createDT <- function(shrinkage.method){
  # create basis with given shrinkage method
  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,shrinkage.method)
  basis.mat.emp[which(is.na(basis.mat.emp))]<-0
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  bb.mat.emp<-create_ds_matrix(jia.DT,shrink.DT,shrinkage.method)
  bb.mat.emp[which(is.na(bb.mat.emp))]<-0
  pred.emp <- predict(pc.emp,newdata=bb.mat.emp)
  emp<-rbind(pc.emp$x,pred.emp)
  emp<-g(emp)[,is.biobank:=FALSE]
  emp[grep("^bb",trait),is.biobank:=TRUE]
  emp[,is.imb.ai:=compare %in% names(ml)]
  emp[,glabel:=trait]
  emp[is.imb.ai==FALSE,glabel:='']
  emp[trait=='control',glabel:='CONTROL']
  emp[trait=='bb_colitis',glabel:='bb_colitis_not_crohns_or_ulcerative_colitis']
  setkey(emp,trait)
  merged <- mani[emp]
  merged[trait=='control',basis.var:=0]
  merged[,PC.ci:=sqrt(basis.var)*1.96]
  list(plot=merged,p1.var=signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3),p2.var=signif(summary(pc.emp)[['importance']][2,]["PC2"]*100,digits=3),e.val=summary(pc.emp)[['importance']][2,])
}



distScores <- function(obj){
  td <- melt(obj$plot,id.vars=c('trait','PC.ci'),measure.vars=paste0('PC',1:11))
  setkeyv(td,c('variable'))
  cont<-td[trait=='control',.(variable,covalue=value)]
  setkeyv(cont,c('variable'))
  merge <- cont[td]
  varv<-data.table(variable=names(obj$e.val),var=obj$e.val)
  setkey(varv,variable)
  merge <- merge[varv]
  merge[,c('dist','wdist','wdists','Zw'):=list(value-covalue,(value-covalue) * var,(value-covalue)^2 * var,((value-covalue)/(PC.ci/1.96)))]
  return(merge)
}

obj<-createDT(SHRINKAGE_METHOD)

dist <- distScores(obj)

pjia<-melt(dist,id.vars=c('trait','variable'),measure.vars='dist') %>% dcast(.,trait ~ variable)

mat <- pjia[,paste0('PC',1:11),with=FALSE] %>% as.matrix
rownames(mat)<-pjia$trait
#pheatmap(mat.old,cluster_cols=FALSE)
plot(hclust(dist(mat)))
plot(hclust(dist(mat[rownames(mat) %in% jia_traits,])))

## non biplot

DT<-melt(obj$plot,id.vars=c('trait','PC.ci'),measure.vars=paste0('PC',1:11))
mins <- DT[,.SD[which.min(value),],by='variable']
maxs <- DT[,.SD[which.max(value),],by='variable']


DT[value %in% c(mins$value,maxs$value),label:=trait]
DT[is.na(label),label:='']
DT[grep('jia',trait),cat:=trait]

library(cowplot)
library(ggrepel)

cDT<-DT[trait=='control',.(variable,control.v=value)]
setkey(DT,'variable')
setkey(cDT,'variable')
DT<-cDT[DT]
DT[,value:=value-control.v]

ggplot(DT,aes(x=variable,y=value,group=trait,label=label,col=cat,alpha=!is.na(cat),ymin=value-PC.ci,ymax=value+PC.ci)) +
geom_line(position = position_dodge(width = 0.9)) + geom_text_repel() + geom_linerange(position = position_dodge(width = 0.9)) +
background_grid(major = "xy", minor = "none")

## can we plot the rotations for PC3 by multiplying the GWAS them ?

rot<-data.table(pid=rownames(pc.emp$rotation),PC3=pc.emp$rotation[,"PC3"])
setkey(rot,pid)
tm<-jia.DT[rot]
tm <- tm[,list(pid=pid,PC3=PC3,p.value=p.value,ppi=wakefield_pp(p.value,maf,unique(n),unique(n1/n))),by=c('trait','ld.block')]
tm[,mscore:=PC3 * ppi]



## as a summary just get the max score for an LD block
forp<-tm[,.SD[which.max(abs(mscore)),],by=c('trait','ld.block')]

ford <- melt(forp,id.vars=c('pid','ld.block','trait'),measure.vars='mscore')
ford <- dcast(ford,ld.block~variable+trait)
rownames(ford)<-ford$ld.block
mat<-as.matrix(ford[,2:10,with=FALSE])
rownames(mat) <- ford$ld.block
library(pheatmap)
pheatmap(mat[apply(mat,1,function(x) max(abs(x))) > 0.01,])

ggplot(forp,aes(x=ld.block,y=mscore)) + geom_point() + facet_wrap(~trait)


if(FALSE){
  basis.mat.est <- create_ds_matrix(basis.DT,shrink.DT,'est')
  ## need to add control where beta is zero
  basis.mat.est<-rbind(basis.mat.est,control=rep(0,ncol(basis.mat.est)))
  pc.est <- prcomp(basis.mat.est,center=TRUE,scale=FALSE)
  jia.mat.est<-create_ds_matrix(jia.DT,shrink.DT,'est')
  pred.emp <- predict(pc.est,newdata=jia.mat.est)
  est<-rbind(pc.emp$x,pred.emp)
  est<-g(est)
  ggplot(est,aes(x=PC1,y=PC2,color=compare,label=trait)) + geom_point() + geom_text() + theme_bw() + ggtitle('Estimate MAF SE shrinkage')
}
