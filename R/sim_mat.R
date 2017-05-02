library(data.table)
library(snpStats)
## this is code to simulate OR under a synthetic null.
(load("/home/ob219/scratch/as_basis/tmp/all_or_shared_with_af.RData"))
## loads into final.t
out.dir<-'/home/ob219/scratch/as_basis/support/1KGenome_snpStats'
## get a uniq list of positions
setkey(final.t,id)
lu<-unique(final.t)[,c('chr','position'),with=FALSE]
## get a list of vcf files
vcf.dir<-'/home/ob219/scratch/DATA/1kgenome/VCF/EUR/by.chr.phase3/'
vcf<-list.files(path=vcf.dir,pattern="ALL.*gz$",full.names=TRUE)
names(vcf)<-sub(".*chr([^\\.]+).*","\\1",vcf)

## use vcftools to select 300K snps we want. This might be super slow.

##get chr snpMatrix object

getSM<-function(DT,vcf.file){
	tmp.file<-tempfile(pattern = "vcf_tmp", tmpdir = tempdir())
	write.table(DT,file=tmp.file,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	## use tabix to get the VCF header
	tabix_cmd<-sprintf("/home/ob219/bin/htslib/tabix -H %s",vcf.file) 
	message(tabix_cmd)
	my.pipe<-pipe(tabix_cmd)
	header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
	close(my.pipe)
	cnames<-unlist(strsplit(header,"\t"))
	vcftools_cmd<-sprintf("/home/ob219/bin/vcftools --gzvcf %s --positions %s --recode -c | grep -v '^#'",vcf.file,tmp.file)
	message(vcftools_cmd)
	tmp<-fread(vcftools_cmd)
	setnames(tmp,cnames)
	unlink(tmp.file)
	gt<-tmp[,10:ncol(tmp),with=FALSE]
	if(nrow(gt)==0)
	  return(NA)
	info<-tmp[,1:9,with=FALSE]
	message("Creating snpMatrix obj")
	sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
	sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
	sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
	## set anything else to a missing value
	sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
	colnames(sm)<-make.names(info$ID,unique=TRUE)
	rownames(sm)<-colnames(gt)
	sm<-new("SnpMatrix", sm)
	return(list(map=info,gt=sm))
}
lapply(split(lu,lu$chr),function(x){
	ofile<-file.path(out.dir,paste0(unique(x$chr),'.RData'))
	vcf.file<-vcf[unique(x$chr)]
	sm<-getSM(x,vcf.file)
	save(sm,file=ofile)
	message(sprintf("Saved %s",ofile))
})

## next add LD assignments
ldBlockGR<-function(file){
        v<-scan(file,"character")
        tmp<-strsplit(gsub("[:-]+",":",v),":")
        tmp<-rbindlist(lapply(tmp,function(x) data.table(chr=x[1],start=x[2],end=x[3])))
        tmp<-tmp[order(as.numeric(tmp$chr),as.numeric(tmp$start)),]
        with(tmp,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
}

ld.gr<-ldBlockGR('/scratch/ob219/as_basis/support/all.1cM.bed')
fs<-list.files(path=out.dir,pattern="*.RData",full.names=TRUE)
lapply(fs,function(f){
	message(sprintf("Processing %s",f))
	load(f)
	## loads into sm object
	m<-sm$map
	setnames(m,'#CHROM','CHROM')
	snp.gr<-with(m,GRanges(seqnames=Rle(CHROM),ranges=IRanges(start=POS,end=POS),id=1:nrow(m)))
	ol<-as.matrix(findOverlaps(snp.gr,ld.gr))
	m$LDBLOCK<-integer(length=nrow(m))
	m[ol[,1],]$LDBLOCK<-ol[,2]
	sm$map<-m
	save(sm,file=f)
})

## next compute a bunch of z scores based on the mvn

library(Matrix)
library(corpcor)
library(mvtnorm)

mvs.perm<-function(sigma,n=1000){
	if(!is.matrix(sigma))
		stop("sigma parameter is not a matrix")		
	if(!is.positive.definite(sigma,,method="chol"))
		stop("sigma is not positive definite")
	## in original paper method="chol" was not defined so I assume used eigen default
	## this is slower than the choleski decomp ! Perhaps we should contact the author ?
	rd<-rmvnorm(n,mean=rep(0,ncol(sigma)),sigma=sigma,method="chol")
	t(rd)
}

attempt.pos.def<-function(mat,diag.val=1.0001){
  print(paste("diag.val",diag.val))
	if(!is(mat,"Matrix"))
		stop("mat is not a Matrix!")
	if(diag.val >= 1.1){
	  print("Matrix is not positive definite. Finding closest approximation..")
		diag(mat)<-1
		return(as(make.positive.definite(mat),"Matrix"))
	}
	diag(mat)<-diag.val
	if(is.positive.definite(mat,,method="chol")==FALSE){
	  new.diag<-signif(1+((diag.val-trunc(diag.val))*10))
		mat<-attempt.pos.def(mat,new.diag)
	}else{
		return(mat)
	}
}

mvs.sigma.r2<-function(r2){
	diag(r2)<-1
	if(!is.positive.definite(r2,,method="chol")){
		#this recurses through various values of diag if we exceed 1 then
		#we compute the closest matrix that is positive definite.
		r2<-attempt.pos.def(r2)
	}
	r2
}

load(f)

ld.idx<-split(1:nrow(sm$map),sm$map$LDBLOCK)
lapply(ld.idx,function(i){
	idx<-ld.idx[[1]]
})
gt<-sm$gt[,idx]
colnames(gt)<-sm$map$POS
tld<-ld(gt,gt,stats="R.squared")
sigma<-as.matrix(mvs.sigma.r2(Matrix(tld)))
perms<-mvs.perm(sigma,n=100)
