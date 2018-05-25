library(cupcake)
library(parallel)
support.dir<-'/scratch/ob219/as_basis/support_tab'
# reference allele frequencies
ref_af_file<-file.path(support.dir,'as_basis_snp_support.tab')
m_file<-file.path(support.dir,'as_manifest_december.tab')
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/scratch/as_basis/gwas_stats/input_files'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,'emp')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)


#input.dir <- '/scratch/ob219/as_basis/jia_ind_analysis//ind_proj_2_0.01_2500_1e+06/split/'

project_individual<-function(input.dir){
  in.files <- list.files(path=input.dir,pattern='*.RDS',full.names=TRUE)
  all.preds <- mclapply(in.files,function(f){
    message(f)
    obj <- readRDS(f)
    snps <- obj$snps
    or <- data.table(exp(obj$proj_lor))
    or[,snp:=snps$pid]
    samples <- obj$samples
    setnames(or,c(as.character(samples$ID_1),'pid'))
    DT <- melt(or,id.vars='pid')
    setnames(DT,c('pid','sample','or'))
    setcolorder(DT,c('sample','pid','or'))
    setkey(DT,pid)
    setnames(DT,'sample','trait')
    mat.emp <- create_ds_matrix(DT,shrink.DT,'emp')
    if(!identical(colnames(mat.emp),colnames(basis.mat.emp)))
      stop("Something wrong basis and projection matrix don't match")
      list(sample=obj$samples,pred=predict(pc.emp,newdata=mat.emp))
    },mc.cores=8)


    all.proj<-data.table(melt(do.call('rbind',lapply(all.preds,'[[','pred'))))
    setnames(all.proj,c('sample','pc','value'))
    all.samples<-rbindlist(lapply(all.preds,'[[','sample'))
    setkey(all.proj,sample)
    setkey(all.samples,ID_1)
    all.proj <- all.proj[all.samples[,.(ID_1,alt_ilar_code)]]
    setnames(all.proj,'alt_ilar_code','trait')
    setcolorder(all.proj,c('trait','pc','value','sample'))


    jia_traits<-fread(m_file)[grep('jia_',trait),]$trait
    #bb_traits<-bb_traits[bb_traits != 'jia_cc']
    jia.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir,jia_traits)
    jia.mat.emp <- create_ds_matrix(jia.DT,shrink.DT,'emp')
    jia.pred <- data.table(melt(predict(pc.emp,newdata=jia.mat.emp)))
    setnames(jia.pred,c('trait','pc','value'))
    jia.pred[,trait:=gsub("jia\\_","",trait)]
    comb<-rbind(jia.pred,all.proj,fill=TRUE)
    comb[,class:=ifelse(is.na(sample),'summary','individual')]

    ## really want to compare control
    control <- data.table(melt(pc.emp$x))[Var1=='control',]
    setnames(control,c('trait','pc','value'))
    control[,class:='control']
    control<-rbindlist(lapply(unique(comb$trait),function(tra){
      tmp<-copy(control)
      tmp[,trait:=tra]
    }))

    final <- rbind(control,comb,fill=TRUE)
    final[is.na(sample),sample:=paste(trait,class,sep='_')]
    final
}


gain_0 <- project_individual('/scratch/ob219/as_basis/jia_ind_analysis//ind_proj_2_0.01_2500_1e+06/split/')
gain_1 <- project_individual('/scratch/ob219/as_basis/jia_ind_analysis//ind_proj_3_0.01_2500_1e+06/split/')
gain_2 <- project_individual('/scratch/ob219/as_basis/jia_ind_analysis//ind_proj_10_0.05_2500_1e+06/split/')

all<-list(gain_0,gain_1,gain_2)

library(ggplot2)
library(cowplot)
 ppl <- ggplot(gain_0[!trait %in% c('cc','missing','UnA'),],aes(x=pc,y=as.numeric(value),alpha=class!='individual',colour=class,group=sample)) +
 geom_point() +  geom_line() + facet_grid(trait~.) + scale_alpha_discrete(guide=FALSE) + xlab("Principle Component") +
 scale_color_manual(name="log(OR) Source",label=c('Control (log(OR)=1)','Individual','Summary'),values=c(control='firebrick',individual='seagreen4',summary='dodgerblue')) +
 theme(legend.position='bottom') +  background_grid(major = "xy", minor = "none") + ylab("Principle component loading")
ppl
save_plot("~/tmp/ind_jia_scree_plot.pdf",ppl,base_height=7)

## show what happens when we turn up the gain !

pc2.plot <- rbindlist(lapply(seq_along(all),function(i){
  d<-all[[i]]
  tmp<-d[trait=='ERA' & pc=='PC2' & class=='individual',]
  tmp[,param:=LETTERS[i]]
  tmp
}))

# get summary and control lines
control.intercept <- gain_0[class=='control' & pc=='PC2' & trait=='ERA',]$value
era.summary.intercept <- gain_0[class=='summary' & pc=='PC2' & trait=='ERA',]$value

ppr<-ggplot(pc2.plot,aes(x=param,y=value,group=param)) + geom_boxplot() + geom_jitter(color='seagreen4',alpha=0.2) + xlab("Posterior log(OR) Constraint") + ylab("ERA PC2 Loading") +
geom_hline(yintercept=control.intercept,size=1,color='firebrick') + geom_hline(yintercept=era.summary.intercept,color='dodgerblue',size=1) +
background_grid(major = "xy", minor = "none") + scale_x_discrete(labels=c("A" = "P(OR>2) = 0.01", "B" = "P(OR>3) = 0.01","C" = "P(OR>10) = 0.05"))
#save_plot("~/tmp/ind_era_boxplot.pdf",ppr,base_height=7)

## if we compare ERA projections of PC2 with projections for other subtypes is that significant ?
## for each parameter setting
all.t <- rbindlist(lapply(seq_along(all),function(i){
  d <- all[[i]]
  ft <- d[class=='individual' & pc=='PC2' & !trait %in% c('cc','missing','UnA'),]
  test.jia.traits <- c('EO','ERA','PO','PsA','RFneg','RFpos','sys')
  rbindlist(lapply(test.jia.traits,function(jt){
    tt <- t.test(ft[trait==jt,]$value,mu=control.intercept)
    data.table(trait=jt,p.val=tt$p.value,parameter=LETTERS[i])
  }))
}))

## do a tileplot type heatmap

ppq <- ggplot(all.t,aes(x=parameter,y=trait,fill=-log(p.val),label=signif(p.val,digits=2))) + geom_tile() + geom_text(color='black') +
scale_x_discrete(labels=c("A" = "P(OR>2) = 0.01", "B" = "P(OR>3) = 0.01","C" = "P(OR>10) = 0.05")) + xlab("Posterior log(OR) Constraint") +
ylab("JIA Subtype") +  scale_fill_gradient(low = "grey", high = "white",guide=FALSE) + theme(axis.line=element_blank())

ppal <- plot_grid(ppr,ppq,nrow=2,labels = c("A", "B"))
save_plot("~/tmp/ind_era_boxplot.pdf",ppal,base_height=8)

ft <- final[class=='individual' & pc=='PC2' & !trait %in% c('cc','missing','UnA'),]
t.test(ft[trait=='ERA',]$value,mu=final[class=='control' & sample=='ERA_control' & pc=='PC2',]$value)

test.jia.traits <- c('EO','ERA','PO','PsA','RFneg','RFpos','sys')
rbindlist(lapply(test.jia.traits,function(jt){
  tt <- t.test(ft[trait==jt,]$value,mu=final[class=='control' & sample=='ERA_control' & pc=='PC2',]$value)
  data.table(trait=jt,p.val=tt$p.value)
}))
