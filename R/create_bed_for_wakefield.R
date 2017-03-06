library(data.table)
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

logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

## compute variance shrinkage for quantitative trait study
Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

## compute variance shrinkage for case control study
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(p,f, N, s,pi_i) {
  sd.prior <- 0.2
  V <- Var.data.cc(f, N, s)
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  sBF <- logsum(lABF + log(pi_i))
  exp(lABF + log(pi_i))/(exp(sBF) + 1)
  #ret <- data.frame(V,z,r,lABF,ppi)
  #if(!is.null(suffix))
  #  colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  #return(ret)
}

ld.gr<-ldBlockGR('/scratch/ob219/as_basis/support/all.1cM.bed')
snp.gr<-with(final.t,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,end=position),id=1:nrow(final.t)))

ol<-as.matrix(findOverlaps(snp.gr,ld.gr))
final.t$ld.block<-0
final.t[ol[,1],]$ld.block<-ol[,2]

by.disease<-split(final.t,final.t$disease)

add.pp<-function(DT,pi_i=1e-4,total,prop){
	pp<-with(DT,approx.bf.p(p.val,maf,total,prop,pi_i))
	DT$pp<-pp
	return(DT)
}

final.t$maf<-final.t$risk.allele.freq
final.t[maj_idx,]$maf<-1-final.t[maj_idx,]$maf

## for each disease 
all.pp<-lapply(names(by.disease),function(n){
	message(paste("Processing",n))
	ss.idx<-which(ss$disease==n)
	cases<-ss[ss.idx,]$cases
	controls<-ss[ss.idx,]$controls
	total<-cases+controls
	prop.case<-signif(cases/total,digits=3)
	tmp<-by.disease[[n]]
	by.ld<-split(tmp,tmp$ld.block)
	pp<-lapply(by.ld,add.pp,total=total,prop=prop.case)
	rbindlist(pp)
})

final.t<-rbindlist(all.pp)
final.t$por<-log(final.t$or) * final.t$pp


## f is defined as the minor allele frequency
maj_idx<-which(final.t$risk.allele.freq>0.5)
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
save(no.pp,file="/home/ob219/scratch/as_basis/tmp/no_p_pp_matrix.RData")
#save(final,file="/home/ob219/scratch/as_basis/tmp/all_or_shared.RData")


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
