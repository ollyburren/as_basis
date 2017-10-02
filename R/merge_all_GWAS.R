library(data.table)
library(GenomicRanges)
## processed summary stats
## FUNCTIONS
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

process.dir<-"/scratch/ob219/as_basis/gwas_stats/processed/"
allele_freq_file <- '/scratch/ob219/as_basis/tmp/bb_all_or_shared_with_af.RData'
final_output_file <- '/home/ob219/scratch/as_basis/merged_data/with_biobank_traits.RData'

if(!file.exists(allele_freq_file)){
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
  if(length(dup.in.cohort) > 0)
	 shared.or<-shared.or[-which(shared.or$id %in% dup.in.cohort),]
	final.keep<-unique(shared.or$id)

	final<-all.id[which(all.id$id %in% final.keep),]
	## prepare final stub so we switch alleles if required
	stub<-final[!duplicated(final$id),c('id','risk.allele','other.allele'),with=FALSE]
	setkey(stub,'id')
	setkey(final,'id')
	setnames(stub,c('id','a1','a2'))
	#merge back with final
	final<-final[stub]
	flip.idx<-with(final,which(risk.allele == a2 & other.allele ==a1))
	final[flip.idx,]$or<-signif(1/final[flip.idx,]$or,digits=3)
	## fix alleles
	final[flip.idx,]$risk.allele<-final[flip.idx,]$a1
	final[flip.idx,]$other.allele<-final[flip.idx,]$a2
	non.match.idx<-which(final$risk.allele != final$a1)
	## are these because they are complement
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
  del.idx<-which(final$id %in% unique(final[which(final$risk.allele != final$a1),]$id))
  if(length(del.idx>0))
	 final<-final[-del.idx,]
  final$position<-as.integer(final$position)
	save(final,file="/scratch/ob219/as_basis/tmp/combined_raw_as_basis.RData")

	## final is a long thin list of OR here is a function to create (based on p-val threshold)
	## a matrix of log(OR) - so they are symmetrical. Also adds a control column  where
	## all OR =1 or log(or)=0

	## add the allele frequencies
	#load("/home/ob219/scratch/as_basis/tmp/all_or_shared.RData")
	load("/home/ob219/scratch/as_basis/1KG_support/all_EUR_support.RData")
  ## grab the stub - this should be the unique list of SNPs for which we want to get AF
  stub[,c('chr','position'):=tstrsplit(stub$id,':')]
  stub$position<-as.integer(stub$position)
  setnames(stub,paste('stub',names(stub),sep='_'))
  setkeyv(stub,c('stub_chr','stub_position'))


	## get just the gt that we need
	#tmp<-unique(final,by='id')[,c('chr','position','risk.allele','other.allele'),with=FALSE]
  ## think this is different
	setkeyv(all.eur,c('chr','position'))
  support<-merge(all.eur,stub,by.x=c('chr','position'),by.y=c('stub_chr','stub_position'))
  support$stub_a2.AF<-double(length=nrow(support))
  ## direct match
  ok.idx<-with(support,which(a2==stub_a2 & a1==stub_a1))
  support[ok.idx,]$stub_a2.AF<-support[ok.idx,]$a2.f
  ## flip
  flip.idx<-with(support,which(a1==stub_a2 & a2==stub_a1))
  support[flip.idx,]$stub_a2.AF<-1-support[flip.idx,]$a2.f
  ## revcom normal
  ncomp.idx<-with(support,which(comp(stub_a1) == a1 & comp(stub_a2) ==a2 & stub_a2.AF==0))
  support[ncomp.idx,]$stub_a2.AF<-support[ncomp.idx,]$a2.f
  ## revomp flipped
  fcomp.idx<-with(support,which(comp(stub_a2) == a1 & comp(stub_a1) ==a2 & stub_a2.AF==0))
  support[fcomp.idx,]$stub_a2.AF<-1-support[fcomp.idx,]$a2.f
  wrong.idx<-setdiff(1:nrow(support),c(ok.idx,flip.idx,ncomp.idx,fcomp.idx))
  if(length(wrong.idx)>0)
    stop("Problem getting reference allele freq")
  ## if we don't get the above means that we agree on alleles and therefore only need to store
  #MAF and the allele with which it goes
	csupport<-unique(support[,c('chr','position','name','stub_a1','stub_a2','stub_a2.AF'),with=FALSE])
  csupport[,minor_allele:=stub_a2]
  csupport[csupport$stub_a2.AF>0.5,minor_allele:=stub_a1]
  csupport[,maf:=stub_a2.AF]
  csupport[csupport$stub_a2.AF>0.5,maf:=1-stub_a2.AF]
  csupport<-csupport[,c('chr','position','name','minor_allele','maf'),with=FALSE]
	save(csupport,file='/scratch/ob219/as_basis/tmp/shared_support_file_with_AF.RData')
	## add these back to final
	final.t<-merge(final,csupport,by.x=c('chr','position'),by.y=c('chr','position'))
  ## get the minor allele
	#final.t$maf<-final.t$risk.allele.freq
	#final.t[final.t$maf>0.5,]$maf<-1-final.t[final.t$maf>0.5,]$maf
	save(final.t,file=allele_freq_file)
}else{
	## loads into final.t
	load(allele_freq_file)
	## compute maf
}


## need to update this
ss<-fread("/home/ob219/scratch/as_basis/gwas_stats/sample_counts_bb.csv")
ss$prop<-0
#by.disease<-split(final.t,final.t$disease)

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
#final.t$ld.block<-0
final.t$ld.block<-integer(length=nrow(final.t))
final.t[ol[,1],]$ld.block<-ol[,2]

by.disease<-split(final.t,final.t$disease)

add.pp<-function(DT,pi_i=1e-4,total,prop,type){
	pp<-with(DT,approx.bf.p(p.val,maf,total,prop,pi_i,type))
	DT$pp<-pp
	return(DT)
}

#maj_idx<-which(final.t$risk.allele.freq>0.5)
#final.t$maf<-final.t$risk.allele.freq
#final.t[maj_idx,]$maf<-1-final.t[maj_idx,]$maf
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
save(final.t,file=final_output_file)
