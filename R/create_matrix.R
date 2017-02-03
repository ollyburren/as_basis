library(ggbiplot)
library(data.table)
process.dir<-"/scratch/ob219/as_basis/gwas_stats/processed/"
adf<-list.files(path=process.dir,pattern='*.RData',full.names=TRUE)

ad<-lapply(seq_along(adf),function(i){
	message(paste("Processing",basename(adf[i])))
	tmp<-get(load(adf[i]))
	tmp$disease<-sub("\\.RData","",basename(adf[i]))
	tmp
})

## work out which SNPs are found in all datasets


all.id<-rbindlist(lapply(ad,function(d){
	d$id<-paste(d$chr,d$position,sep=':')
	d
}))

all.f<-all.id[,c('id','or','disease','p.val'),with=FALSE]
setkey(all.f,id)
## which SNPs are found across all cohorts
id.count<-table(all.id$id)
keep.id<-names(id.count[id.count==length(unique(all.f$disease))])
keep.idx<-which(all.f$id %in% keep.id)
shared.or<-all.f[keep.idx,]
## there are some duplicates here that are bumping up the numbers
foo<-split(shared.or$id,shared.or$disease)
dup.in.cohort<-unique(do.call('c',lapply(foo,function(id) id[duplicated(id)])))
shared.or<-shared.or[-which(shared.or$id %in% dup.in.cohort),]
final.keep<-unique(shared.or$id)

final<-all.id[which(all.id$id %in% final.keep),]
## prepare final stub so we switch alleles if required
stub<-final[!duplicated(final$id),c('id','risk.allele','other.allele'),with=FALSE]
setkey(stub,'id')
setkey(final,'id')
setnames(stub,c('id','a1','a2'))
final<-final[stub]
flip.idx<-with(final,which(risk.allele == a2 & other.allele ==a1))
final[flip.idx,]$or<-signif(1/final[flip.idx,]$or,digits=3)
## fix alleles
final[flip.idx,]$risk.allele<-final[flip.idx,]$a1
final[flip.idx,]$other.allele<-final[flip.idx,]$a2
non.match.idx<-which(final$risk.allele != final$a1)
## are these because they are complement 
comp<-function(cv){
  tmp<-cv
  a<-c('A','G','C','T')
  b<-c('T','C','G','A')
  idx<-split(1:length(cv),cv)
  for(i in seq_along(a)){
    cv[idx[[a[i]]]]<-b[i]
  }
  cv	
}

## how many have alleles where it is impossible to work out where rev comp e.g.
#A,T T,A G,C C,G -- there is one and it's already taken care of
#subset(final,risk.allele=='C' & other.allele=='G')
flip.comp.idx<-with(final,which(comp(risk.allele) == a2 & comp(other.allele) ==a1))
##change 
actual.flip.comp.idx<-intersect(non.match.idx,flip.comp.idx)
final[actual.flip.comp.idx,]$or<-signif(1/final[actual.flip.comp.idx,]$or,digits=3)
## fix alleles
final[actual.flip.comp.idx,]$risk.allele<-final[actual.flip.comp.idx,]$a1
final[actual.flip.comp.idx,]$other.allele<-final[actual.flip.comp.idx,]$a2
## these are just where they have used complement alleles rather than a switch
rev.comp.idx<-with(final,which(comp(risk.allele) == a1 & comp(other.allele) ==a2))
final[rev.comp.idx,]$risk.allele<-final[rev.comp.idx,]$a1
final[rev.comp.idx,]$other.allele<-final[rev.comp.idx,]$a2
## left overs are odd indel and should be deleted
final<-final[-which(final$id %in% unique(final[which(final$risk.allele != final$a1),]$id)),]

## final is a long thin list of OR here is a function to create (based on p-val threshold)
## a matrix of log(OR) - so they are symmetrical. Also adds a control column  where 
## all OR =1 or log(or)=0

createORMatrix<-function(DT,p.val.thresh=1){
   tmp<-unique(DT[DT$p.val<p.val.thresh,]$id)
   tmp<-DT[which(DT$id %in% tmp),]
   tmp$lor<-log(tmp$or)
   tmp<-melt(tmp,id.vars=c('id','disease'),measure.vars = 'lor') 
   ret<-dcast(tmp,disease~id)
   diseases<-ret$disease
   ret<-as.data.frame(ret[,2:ncol(ret),with=FALSE])
   rownames(ret)<-diseases
   fret<-rbind(ret,rep(0,ncol(ret)))
   rownames(fret)<-c(diseases,'control')
   fret
}

no.p<-createORMatrix(final)
save(no.p,file="/home/ob219/scratch/as_basis/tmp/no_p_lor_matrix.RData")
save(final,file="/home/ob219/scratch/as_basis/tmp/all_or_shared.RData")


### p.vals cut offs 
#p.val<-c(1,0.1,1e-2,1e-3,1e-4)
#
### what about crazy or > 10
#
#
#function(p){
#   tmp<-unique(final[final$p.val<p,]$id)
#   tmp<-final[which(final$id %in% tmp),]
#   tmp$lor<-log(tmp$or)
#   tmp<-melt(tmp,id.vars=c('id','disease'),measure.vars = 'lor') 
#   ret<-dcast(tmp,disease~id)
#   diseases<-ret$disease
#   ret<-as.data.frame(ret[,2:ncol(ret),with=FALSE])
#   rownames(ret)<-diseases
#   ret.pca <- prcomp(ret,center=TRUE,scale=TRUE)
#   g <- ggbiplot(ret.pca, choices = 2:3, obs.scale = 1, var.scale = 1,labels=diseases,
#               ellipse = TRUE, circle = TRUE,var.axes=FALSE)
#   g <- g + coord_cartesian(xlim=c(-10,10),ylim=c(-10,10))
#   
#}
#
#save(ret,file="/home/ob219/scratch/as_basis/tmp/lor_matrix.RData")
