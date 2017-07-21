library(data.table)
library(GenomicRanges)
process.dir<-"/scratch/ob219/as_basis/gwas_stats/processed/"

if(!file.exists("/scratch/ob219/as_basis/tmp/all_or_shared_with_af.RData")){
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

## add the allele frequencies
#load("/home/ob219/scratch/as_basis/tmp/all_or_shared.RData")
load("/home/ob219/scratch/as_basis/1KG_support/all_EUR_support.RData")
## get just the gt that we need
tmp<-unique(final)[,c('chr','position','risk.allele','other.allele'),with=FALSE]
setkey(all.eur,chr,position)
setkey(tmp,chr,position)
support<-merge(all.eur,tmp,by.x=c('chr','position'),by.y=c('chr','position'))
## next check alleles
support$risk.allele.freq<-double(length=nrow(support))
ok.idx<-with(support,which(a2==risk.allele & a1==other.allele))
support[ok.idx,]$risk.allele.freq<-support[ok.idx,]$a2.f
flip.idx<-with(support,which(a1==risk.allele & a2==other.allele))
support[flip.idx,]$risk.allele.freq<-1-support[flip.idx,]$a2.f
csupport<-support[,c('chr','position','name','risk.allele.freq'),with=FALSE]
save(csupport,file='/scratch/ob219/as_basis/tmp/shared_support_file_with_AF.RData')
wrong.idx<-setdiff(1:nrow(support),c(ok.idx,flip.idx))
## add these back to final
final.t<-merge(final,csupport,by.x=c('chr','position'),by.y=c('chr','position'))
save(final.t,file='/scratch/ob219/as_basis/tmp/all_or_shared_with_af.RData')
}else{
	## loads into final.t
	load("/scratch/ob219/as_basis/tmp/all_or_shared_with_af.RData")
	## compute maf
	final.t$maf<-final.t$risk.allele.freq
	final.t[final.t$maf>0.5,]$maf<-1-final.t[final.t$maf>0.5,]$maf
}

ss<-fread("/home/ob219/scratch/as_basis/gwas_stats/sample_counts.csv")
ss$prop<-0
by.disease<-split(final.t,final.t$disease)

ldBlockGR<-function(file){
	v<-scan(file,"character")
	tmp<-strsplit(gsub("[:-]+",":",v),":")
	tmp<-rbindlist(lapply(tmp,function(x) data.table(chr=x[1],start=x[2],end=x[3])))
	tmp<-tmp[order(as.numeric(tmp$chr),as.numeric(tmp$start)),]
	with(tmp,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
}

## contains code for computing wakefields aBF and thus posterior probabilities
source("~/git/as_basis/R/wakefield.R")


ld.gr<-ldBlockGR('/scratch/ob219/as_basis/support/all.1cM.bed')
snp.gr<-with(final.t,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,end=position),id=1:nrow(final.t)))

ol<-as.matrix(findOverlaps(snp.gr,ld.gr))
final.t$ld.block<-0
final.t[ol[,1],]$ld.block<-ol[,2]


add.pp<-function(DT,pi_i=1e-4,total,prop,type){
	pp<-with(DT,approx.bf.p(p.val,maf,total,prop,pi_i,type))
	DT$pp<-pp
	return(DT)
}

maj_idx<-which(final.t$risk.allele.freq>0.5)
final.t$maf<-final.t$risk.allele.freq
final.t[maj_idx,]$maf<-1-final.t[maj_idx,]$maf
by.disease<-split(final.t,final.t$disease)
## for each disease 
all.pp<-lapply(names(by.disease),function(n){
	message(paste("Processing",n))
	ss.idx<-which(ss$disease==n)
	cases<-ss[ss.idx,]$cases
	controls<-ss[ss.idx,]$controls
	total<-cases+controls
	prop.case<-signif(cases/total,digits=3)
	t<-'CC'
	if(prop.case==1)
	  t<-'QUANT'
	tmp<-by.disease[[n]]
	by.ld<-split(tmp,tmp$ld.block)
	pp<-lapply(by.ld,add.pp,total=total,prop=prop.case,type=t)
	rbindlist(pp)
})

final.t<-rbindlist(all.pp)
final.t$lor<-log(final.t$or)
final.t$por<-log(final.t$or) * final.t$pp

final.t$Z<-sign(final.t$lor) * qnorm(final.t$p.val/2,lower.tail=FALSE)
final.t$pZ<-final.t$Z * final.t$pp

## final thing to try is adjusting Z so that includes variance due to sample size

by.disease<-split(final.t,final.t$disease)
all.adjZ<-lapply(names(by.disease),function(n){
        message(paste("Processing",n))
        ss.idx<-which(ss$disease==n)
        cases<-ss[ss.idx,]$cases
        controls<-ss[ss.idx,]$controls
        total<-cases+controls
        tmp<-by.disease[[n]]
	ret<-tmp$Z * sqrt(1/total)
	#if(controls==0){
	#	ret<-tmp$Z * sqrt(1/total)
	#}else{
	#	ret<-tmp$Z * sqrt((1/cases) + (1/controls))
	#}
	return(ret)
})

final.t$Zadj<-do.call('c',all.adjZ)
final.t$pZadj<-final.t$Zadj * final.t$pp

save(final.t,file='/scratch/ob219/as_basis/tmp/final.t.RData')
## f is defined as the minor allele frequency
##split by sample and compute the partial variance 
## this should be done on a study by study basis

createORMatrix<-function(DT,p.val.thresh=1,var='lor'){
   tmp<-unique(DT[DT$p.val<p.val.thresh,]$id)
   tmp<-DT[which(DT$id %in% tmp),]
   tmp$lor<-log(tmp$or)
   tmp<-melt(tmp,id.vars=c('id','disease'),measure.vars = var) 
   ret<-dcast(tmp,disease~id)
   diseases<-ret$disease
   ret<-as.data.frame(ret[,2:ncol(ret),with=FALSE])
   rownames(ret)<-diseases
   fret<-rbind(ret,rep(0,ncol(ret)))
   rownames(fret)<-c(diseases,'control')
   fret
}

no.pp<-createORMatrix(final.t,var='por')
save(no.pp,file="/home/ob219/scratch/as_basis/tmp/no_p_pp_matrix3.RData")
no.lor<-createORMatrix(final.t)
save(no.lor,file="/home/ob219/scratch/as_basis/tmp/no_p_lor_matrix3.RData")
no.Z<-createORMatrix(final.t,var='Z')
save(no.Z,file="/home/ob219/scratch/as_basis/tmp/no_p_Z_matrix.RData")
no.pZ<-createORMatrix(final.t,var='pZ')
save(no.pZ,file="/home/ob219/scratch/as_basis/tmp/no_p_pZ_matrix.RData")
no.Zadj<-createORMatrix(final.t,var='Zadj')
save(no.Zadj,file="/home/ob219/scratch/as_basis/tmp/no_p_Zadj_matrix.RData")
no.pZadj<-createORMatrix(final.t,var='pZadj')
save(no.pZadj,file="/home/ob219/scratch/as_basis/tmp/no_p_pZadj_matrix.RData")


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
