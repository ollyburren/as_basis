library(data.table)
library(snpStats)
library(GenomicRanges)
## this is code to simulate OR under a synthetic null.
#(load("/home/ob219/scratch/as_basis/tmp/all_or_shared_with_af.RData"))
(load("/home/ob219/scratch/as_basis/tmp/final.t.RData"))
## loads into final.t
out.dir<-'/home/ob219/scratch/as_basis/support/simulations/1KGenome_snpStats'
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

## next add the z score of our target GWAS

target.gwas<-'ill.t1d'
tgwas.DT<-subset(final.t,disease==target.gwas)[,c('id','Z'),with=FALSE]
## add in Z scores to files.
of<-list.files(path=out.dir,pattern="*.RData",full.names=TRUE)
for(f in of){
	message(sprintf("Processing %s",f))
	(load(f))
	DT<-data.table(sm$map)
	setnames(DT,'#CHROM','CHROM')
	DT$id<-paste(DT$CHROM,DT$POS,sep=':')
	DT$order<-1:nrow(DT)
	setkey(DT,id)
	setkey(tgwas.DT,id)
	tmp<-tgwas.DT[DT]
	tmp<-tmp[order(tmp$order),]
	tmp$order<-NULL
	sm$map<-tmp
	save(sm,file=f)
}



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

mvs.perm<-function(m,sigma,n=1000){
	## if m is 1 then we just simulate from normal distribution with mean(Z)
	if(length(m)==1){
		return(rnorm(n,mean=m))
	}
	if(!is.matrix(sigma))
		stop("sigma parameter is not a matrix")		
	if(!is.positive.definite(sigma,,method="chol"))
		stop("sigma is not positive definite")
	## in original paper method="chol" was not defined so I assume used eigen default
	## this is slower than the choleski decomp ! Perhaps we should contact the author ?
	rd<-rmvnorm(n,mean=m,sigma=sigma,method="chol")
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
	if(any(is.na(r2)))
		message(sprintf("Found %s where r2 was na",sum(is.na(r2))))
		r2[is.na(r2)]<-0
	return(as(make.positive.definite(r2),"Matrix"))
	if(!is.positive.definite(r2,,method="chol")){
		#this recurses through various values of diag if we exceed 1 then
		#we compute the closest matrix that is positive definite.
	  	print("Matrix is not positive definite. Finding closest approximation..")
		r2<-attempt.pos.def(r2)
	}
	r2
}

## for each chromosome and LD block create 200 simulations
sim.by.chr<-lapply(of,function(f){
	message(sprintf("Processing %s",f))
	sm<-get(load(f))
	ld.idx<-split(1:nrow(sm$map),sm$map$LDBLOCK)
	by.chr<-lapply(seq_along(ld.idx),function(i){
		idx<-ld.idx[[i]]
		message(sprintf("Processing %s",length(idx)))
		if(length(idx)==1){
			message("Only one sample from normal distro")
			## just sample from norm
			return(t(rnorm(200,mean=sm$map[idx,]$Z)))
		}
		gt<-sm$gt[,idx]
		colnames(gt)<-sm$map[idx,]$POS
		tld<-ld(gt,gt,stats="R.squared")
		sigma<-as.matrix(mvs.sigma.r2(Matrix(tld)))
		return(mvs.perm(sm$map[idx,]$Z,sigma,n=200))
	})
	ret<-data.table(do.call('rbind',by.chr))
	ret$id<-sm$map$id
	ret
})

all.chr.sim<-do.call('rbind',sim.by.chr)
## make sure in the same format as pca basis matrix by reshaping the same
library(reshape2)
## long thin
tmp<-melt(all.chr.sim,id.vars=c('id'))
## compute pAdjZ (the posterior prob adjusted Z score ppi * Z * 1/sqrt(N)
ppi<-ill.id[,c('pp','id'),with=FALSE]
tmpy<-merge(tmp,ppi,by.x='id',by.y='id')
## constant variance due to sample size
var_sam<-1/sqrt(3983+3999)
## THIS IS WRONG THE PP needs to be computed from scratch - cannot use tmpy pp
tmpy$nv<-tmpy$value * var_sam * tmpy$pp
tmpy$value<-tmpy$nv
tmpy<-tmpy[,1:3,with=FALSE]
tmpy.sim<-dcast(tmpy,variable~id)
tmpy.sim<-tmpy.sim[,2:ncol(tmpy.sim)]
## split into two sets of simulations
sim<-as.matrix(tmpy.sim[1:100,])
proj<-as.matrix(tmpy.sim[101:200,])



## there are some variants missing - for this analysis it should not matter too much as long as it's not an error.
## instead we need to make sure that these agree with the actual scores we are using for the pca basis.
filt<-subset(final.t,id %in% all.chr.sim$id & disease %in% c('asthma','CD','CEL','eosinophil','JIA_nosys','lymphocyte','MS','myeloid','PBC','PSC','RA','SLE','UC','wbc'))

## create the matrix

createORMatrix<-function(DT,var='lor'){
   DT<-melt(DT,id.vars=c('id','disease'),measure.vars = var) 
   ret<-dcast(DT,disease~id)
   print(class(ret))
   diseases<-ret$disease
   ret<-as.data.frame(ret[,2:ncol(ret)])
   rownames(ret)<-diseases
   fret<-rbind(ret,rep(0,ncol(ret)))
   rownames(fret)<-c(diseases,'control')
   fret
}

nEucledian<-function(simLoad,actLoad,vexp){
    sqrt((actLoad-simLoad)^2 %*% vexp)
}

computeVarWEuc<-function(pc,proj.pc){
    act.pc<-pc$x
    all.pc<-rbind(pc$x,proj.pc)
    ## what is the variance explained
    vexp<-summary(pc)$importance[2,]
    idx<-nrow(all.pc)
    apply(all.pc,1,nEucledian,all.pc[idx,],vexp)["sim"]
}

basis<-createORMatrix(filt,var='pZadj')
save(list=c(basis,proj,sim),file="/home/ob219/scratch/as_basis/tmp/basis_distance.RData")
## next we add a row from sim and compute pca then project from proj
## will take approx an hour to compute serially
lapply(1:100,function(i){
	message(sprintf("Processing %d",i))
	mat<-rbind(basis,sim[i,])
	rownames(mat)[length(rownames(mat))]<-'sim'
	pca<-prcomp(mat,center=TRUE,scale=FALSE)
	proj.pc<-predict(pca,t(proj[i,]))
	computeVarWEuc(pca,proj.pc)
})



## checking that Z * 1/sqrt(n) approximates beta/var_maf

ill.id<-subset(final.t,disease=='ill.t1d')

# a is a1 counts ctl
# b is a2 counts ctl
# c is a1 counts cas
# d is a2 counts cas

## calc a a=n0(1-f)

ca<-function(n0,f){
    n0*(1-f)
}

cb<-function(n0,f){
    n0*f
}

cc<-function(n1,a,b,theta){
    (n1*a)/(a+(b*theta))
}

cd<-function(n1,a,b,theta){
    (n1*b)/(a+(b*theta))
}

calc.pv<-function(n0,n1,f,theta){
    a<-ca(n0,f)
    b<-cb(n0,f)
    c<-cc(n1,a,b,theta)
    d<-cd(n1,a,b,theta)
    recip.sm<-do.call('cbind',lapply(list(a,b,c,d),function(fi) 1/fi))
    return(sqrt(rowSums(recip.sm)))
}

ill.id$vmaf<-calc.pv(3999,3983,ill.id$maf,ill.id$or)

## gamma_hat can be computed two ways log(or)/vmaf or Z * 1/sqrt(n0+n1)

t<-data.table(with(ill.id,cbind(log(or)/vmaf,Z * 1/sqrt(3999+3983))))
library(ggplot2)
ggplot(t[sample(1:nrow(t),1000),],aes(x=V1,y=V2)) + geom_point() + xlab("vmaf") + ylab("Z times vsample")
## roughly the same just different scaling - which means that we should not mix how we calculate. If we want to include QTL then suggest Z * 1/sqrt(samplesize)


