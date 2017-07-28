library(data.table)
library(snpStats)
library(GenomicRanges)


target.gwas<-'aff.t1d'

(load("/home/ob219/scratch/as_basis/merged_data/17_traits.RData"))
## loads into final.t
out.dir<-'/home/ob219/scratch/as_basis/support/simulations/1KGenome_snpStats'
## get a uniq list of positions
setkey(final.t,id)
lu<-unique(final.t)[,c('chr','position'),with=FALSE]
## get a list of vcf files
vcf.dir<-'/home/ob219/scratch/DATA/1kgenome/VCF/EUR/by.chr.phase3/'
vcf<-list.files(path=vcf.dir,pattern="ALL.*gz$",full.names=TRUE)
names(vcf)<-sub(".*chr([^\\.]+).*","\\1",vcf)

out.dir<-file.path(out.dir,target.gwas)
dir.create(out.dir)
## estimate z score from p.value
final.t$Z<-qnorm(0.5 * final.t$p.val, lower.tail = FALSE) * sign(final.t$lor)

## use vcftools to select 300K snps we want. This might be super slow.

##get chr snpMatrix object
of<-list.files(path=out.dir,pattern="*.RData",full.names=TRUE)
if(length(of)==0){
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
	setnames(info,'#CHROM','CHROM')
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

blash<-lapply(split(lu,lu$chr),function(x){
	ofile<-file.path(out.dir,paste0(unique(x$chr),'.RData'))
	vcf.file<-vcf[unique(x$chr)]
	sm<-getSM(x,vcf.file)
	save(sm,file=ofile)
	message(sprintf("Saved %s",ofile))
})


## next add the z score of our target GWAS
tgwas.DT<-subset(final.t,disease==target.gwas)[,c('id','Z','or','maf'),with=FALSE]
tgwas.DT$lor<-with(tgwas.DT,log(or))
## add in Z scores to files.
of<-list.files(path=out.dir,pattern="*.RData",full.names=TRUE)
for(f in of){
	message(sprintf("Processing %s",f))
	(load(f))
	DT<-data.table(sm$map)
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
blash<-lapply(fs,function(f){
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
}else{
	message(sprintf("Already done %s",target.gwas))
}
