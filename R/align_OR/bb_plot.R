library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp'
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_feb/'

## load data

support.dir<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
# reference allele frequencies
ref_af_file<-'/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab'
#ld_file<-file.path(support.dir,'all.1cM.tab')
#m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv'
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/rds/user/ob219/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)

basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

## compute the variances

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


bb_traits<-fread(m_file)[grep('bb_',trait),]$trait
#bb_traits<-c(bb_traits,'breast_cancer')
bb.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,TRUE,bb_traits)
bb.mat.emp<-create_ds_matrix(bb.DT,shrink.DT,SHRINKAGE_METHOD)
bb.mat.emp[which(is.na(bb.mat.emp))]<-0
pred.emp <- predict(pc.emp,newdata=bb.mat.emp)
emp<-rbind(pc.emp$x,pred.emp)
emp <- cbind(as.data.table(emp),trait=rownames(emp))

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
  asthma = 'bb_asthma'
)
merged$compare <- "none"


merged[,glabel:=trait]
merged[,glabel:=gsub("jia\\_","",trait)]
#setnames(merged,'variable','pc')
for(i in seq_along(ml)) {
    merged[trait %in% c(names(ml)[i], ml[i]), compare:=names(ml)[i]]
}
merged[trait=="control",compare:="control"]

plot <- merged[,.(trait,glabel,compare,basis_trait,pc,value,ci)]
pc <- dcast(plot[pc %in% c('PC1','PC2'),],trait+glabel+compare+basis_trait~pc)
ci<-melt(plot[pc %in% c('PC1','PC2'),.(glabel,pc,ci)],id.vars=c('glabel','pc')) %>% dcast(.,glabel~pc)
setnames(ci,c('glabel','ci_PC1','ci_PC2'))
plot <- cbind(pc,ci[,.(ci_PC1,ci_PC2)])
plot[compare=="none",glabel:='']

## want to show some traits event though they have no compare

plot[trait=='bb_hyperthyroidism_thyrotoxicosis',glabel:=trait]
plot[trait=='bb_hypothyroidism_myxoedema',glabel:=trait]
plot[trait=='bb_psoriatic_arthropathy',glabel:=trait]
plot[trait=='bb_colitis',glabel:=trait]
plot[trait=='bb_AS',glabel:=trait]
## biplot

p1.var=signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3)
p2.var=signif(summary(pc.emp)[['«importance']][2,]["PC2"]*100,digits=3)


library(cowplot)
point.size<-1
pp <- ggplot(plot,aes(x=PC1,y=PC2,color=compare,label=glabel,alpha=glabel!='',ymin=PC2-ci_PC2,ymax=PC2+ci_PC2,xmin=PC1-ci_PC1,xmax=PC1+ci_PC2)) +
geom_point(size=point.size) + geom_text_repel()  + scale_color_discrete(guide=FALSE) +
scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + xlab(sprintf("%s (%.1f%%)",'PC1',p1.var)) +
ylab(sprintf("%s (%.1f%%)",'PC2',p2.var)) +  background_grid(major = "xy", minor = "none") +
geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5)
save_plot("~/tmp/bb_talk.pdf",pp,base_aspect_ratio=1,base_width=7,base_height=7)



# ## load ci model
# mani <- fread(m_file)[,.(trait,cases,controls)]
# mani <- mani[,c('total','ncases'):=list(cases+controls,cases)]
# #mani[,l.ss:=log(cases+controls)]
# ## load in model for variance callibration
# callibration <- readRDS("/home/ob219/rds/hpc-work/as_basis/callibration/callibration_model_fixed.RDS")
#
# mani[,basis.var:=exp(predict(callibration,mani))]
# setkey(mani,trait)
#
# ## comparisons to make
# ml<-list(
#   CD = 'bb_CD',
#   CEL = 'bb_CEL',
#   MS = 'bb_MS',
#   RA = 'bb_RA',
#   SLE = 'bb_SLE',
#   T1D = 'bb_T1D',
#   UC = 'bb_UC',
#   PBC = 'PBC',
#   PSC = 'PSC',
#   asthma = 'bb_asthma'
# )
#
#
# ## function for colouring in pairs of diseases
# g <- function(M){
#     M <- cbind(as.data.table(M),trait=rownames(M))
#     M$compare<-"none"
#     for(i in seq_along(ml)) {
#         M[trait %in% c(names(ml)[i], ml[i]), compare:=names(ml)[i]]
#     }
#     M[trait=="control",compare:="control"]
#     M
# }
#
# ## code to project and create plotting data frame with different shrinkage methods
#
# createDT <- function(shrinkage.method){
#   # create basis with given shrinkage method
#   basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,shrinkage.method)
#   basis.mat.emp[which(is.na(basis.mat.emp))]<-0
#   ## need to add control where beta is zero
#   basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
#   pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
#   bb.mat.emp<-create_ds_matrix(bb.DT,shrink.DT,shrinkage.method)
#   bb.mat.emp[which(is.na(bb.mat.emp))]<-0
#   pred.emp <- predict(pc.emp,newdata=bb.mat.emp)
#   emp<-rbind(pc.emp$x,pred.emp)
#   emp<-g(emp)[,is.biobank:=FALSE]
#   emp[grep("^bb",trait),is.biobank:=TRUE]
#   emp[,is.imb.ai:=compare %in% names(ml)]
#   emp[,glabel:=trait]
#   emp[is.imb.ai==FALSE,glabel:='']
#   emp[trait=='control',glabel:='CONTROL']
#   emp[trait=='bb_colitis',glabel:='bb_colitis_not_crohns_or_ulcerative_colitis']
#   setkey(emp,trait)
#   merged <- mani[emp]
#   merged[trait=='control',basis.var:=0]
#   merged[,PC.ci:=sqrt(basis.var)*1.96]
#   list(plot=merged,p1.var=signif(summary(pc.emp)[['importance']][2,]["PC1"]*100,digits=3),p2.var=signif(summary(pc.emp)[['importance']][2,]["PC2"]*100,digits=3),e.val=summary(pc.emp)[['importance']][2,])
# }
#
#
# ## make a plot for old style
#
#
# old.obj <- createDT('emp')
# new.obj <- createDT('ws_emp')

## add in callibration


# point.size<-1
# ppa<-ggplot(old.obj$plot,aes(x=PC1,y=PC2,color=compare,label=glabel,alpha=glabel!='',ymin=PC2-PC.ci,ymax=PC2+PC.ci,xmin=PC1-PC.ci,xmax=PC1+PC.ci)) +
# geom_point(size=point.size) + geom_text_repel()  + scale_color_discrete(guide=FALSE) +
# scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + xlab(sprintf("%s (%.1f%%)",'PC1',old.obj$p1.var)) +
# ylab(sprintf("%s (%.1f%%)",'PC2',old.obj$p2.var)) +  background_grid(major = "xy", minor = "none") +
# geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5)
#
# ppb <- ggplot(old.obj$plot,aes(x=PC1,y=PC2,color=compare,label=trait,,ymin=PC2-PC.ci,ymax=PC2+PC.ci,xmin=PC1-PC.ci,xmax=PC1+PC.ci)) + geom_point(size=point.size) +
# geom_text_repel()  + scale_color_discrete(guide=FALSE)  + coord_cartesian(xlim=c(-0.05,0.05),ylim=c(-0.04,0)) +
# xlab(sprintf("%s (%.1f%%)",'PC1',old.obj$p1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',old.obj$p2.var)) +
# background_grid(major = "xy", minor = "none") +
# geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5)
#
#
#
# ppc<-ggplot(new.obj$plot,aes(x=PC1,y=PC2,color=compare,label=glabel,alpha=glabel!='',,ymin=PC2-PC.ci,ymax=PC2+PC.ci,xmin=PC1-PC.ci,xmax=PC1+PC.ci)) +
# geom_point(size=point.size) + geom_text_repel()  + scale_color_discrete(guide=FALSE) +
# scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) + xlab(sprintf("%s (%.1f%%)",'PC1',new.obj$p1.var)) +
# ylab(sprintf("%s (%.1f%%)",'PC2',new.obj$p2.var)) +  background_grid(major = "xy", minor = "none") +
# geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5)
#
# ppd <- ggplot(new.obj$plot,aes(x=PC1,y=PC2,color=compare,label=trait,,ymin=PC2-PC.ci,ymax=PC2+PC.ci,xmin=PC1-PC.ci,xmax=PC1+PC.ci)) + geom_point(size=point.size) +
# geom_text_repel()  + scale_color_discrete(guide=FALSE)  + coord_cartesian(xlim=c(-0.05,0.05),ylim=c(0,-0.04)) +
# xlab(sprintf("%s (%.1f%%)",'PC1',new.obj$p1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',new.obj$p2.var)) +
# background_grid(major = "xy", minor = "none") +
# geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5)
#
# ppold <-plot_grid(ppa, ppb,labels = c("A", "B"))
# ppnew <-plot_grid(ppc,ppd,labels = c("A", "B"))
# save_plot("~/tmp/old_bb_plot.pdf",ppold,ncol = 2,base_height=5)
# save_plot("~/tmp/new_bb_plot.pdf",ppnew,ncol = 2,base_height=5)

## look at hclust using euc distance


emp<-rbind(pc.emp$x,pred.emp)

controlv <- merged[trait=='control',.(pc,covalue=value)]
td <- merged[,.(trait,pc,value,compare,ci)]
setkey(td,'pc')
setkey(controlv,'pc')
M <- controlv[td]
## add variance explained
vexp<-summary(pc.emp)$importance[2,]
vexp <- data.table(pc=names(vexp),vexp=vexp)
setkey(vexp,'pc')
M <- vexp[M]
M[,dist:=(value-covalue)]
dist<-melt(M,id.vars=c('trait','pc'),measure.vars='dist') %>% dcast(.,trait ~ pc)
foo <- as.matrix(dist[,-1])
rownames(foo)<-dist$trait
hc <- hclust(dist(foo))
cols <- rep("black",length(hc$labels))
cols[grep("^bb_",hc$labels)]<-'olivedrab4'

par(mar=c(0, 4, 4, 2)) # c(bottom, left, top, right)
pdf("~/tmp/bb_tree.pdf",height=7,width=7)
plot(hc, hang=-1,xlab="", sub="",col=cols,cex=0.8,main="")
dev.off()
# as.phylo(hc)
#
# colors = c("red", "dodgerblue", "olivedrab4","firebrick","black","red")
# clus5 = cutree(hc, 6)
# plot(as.phylo(hc), type = "fan", tip.color = colors[clus5])
#      label.offset = 1, cex = 0.7)
#
# plot(hclust(dist(foo)))
#
#
# test<-new.obj$plot
#
# distScores <- function(obj){
#   td <- melt(obj$plot,id.vars=c('trait','PC.ci'),measure.vars=paste0('PC',1:11))
#   setkeyv(td,c('variable'))
#   cont<-td[trait=='control',.(variable,covalue=value)]
#   setkeyv(cont,c('variable'))
#   merge <- cont[td]
#   varv<-data.table(variable=names(obj$e.val),var=obj$e.val)
#   setkey(varv,variable)
#   merge <- merge[varv]
#   merge[,c('dist','wdist','wdists','Zw'):=list(value-covalue,(value-covalue) * var,(value-covalue)^2 * var,((value-covalue)/(PC.ci/1.96)))]
#   return(merge)
# }
#
# old.dist <- distScores(old.obj)
# nd1<-old.dist[,list(maxZ.old=max(Zw*var)),by='trait']
# #nd2<-old.dist[,list(var.exp.old=var,sumZ.old=sum(abs(Zw),na.rm=TRUE)),by=c('variable','trait')]
# new.dist <- distScores(new.obj)
# nd2<-new.dist[,list(maxZ.new=max(Zw*var)),by='trait']
# res<-cbind(nd1,nd2)[,.(trait,maxZ.old,maxZ.new)]
# res <- res[,delta:=maxZ.old-maxZ.new]
#
# pold<-melt(old.dist,id.vars=c('trait','variable'),measure.vars='dist') %>% dcast(.,trait ~ variable)
# pnew<-melt(new.dist,id.vars=c('trait','variable'),measure.vars='dist') %>% dcast(.,trait ~ variable)
#
#
#
#
# library(pheatmap)
#
# mat.old <- pold[,paste0('PC',1:11),with=FALSE] %>% as.matrix
# rownames(mat.old)<-pold$trait
# #pheatmap(mat.old,cluster_cols=FALSE)
# plot(hclust(dist(mat.old)))
#
# mat.new <- pnew[,paste0('PC',1:11),with=FALSE] %>% as.matrix
# rownames(mat.new)<-pnew$trait
# #pheatmap(mat.old,cluster_cols=FALSE)
# plot(hclust(dist(mat.new)))
#
# ## plot for chris for new
#
# ggplot(merge,aes(x=trait,y=Zw)) + geom_bar(stat="identity") + facet_wrap(~variable)
