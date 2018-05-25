## SIMULATION WORKBOOK
#Misc code and notes for t1d_bootstrap simulation

## check the current imputed objects and get a list of available SNPs
## MUST BE DONE ON CSD3

# chromosome 1 failed ;(

# In meantime take a look weighting

library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(Gviz)
library(biomaRt)

args <- list(
    snp_support_file='/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab',
    outfile="~/tmp/testbasis.RDS",
    manifest_file='/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv',
    gwas_data_dir='/rds/user/ob219/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
)

## loaf in basis data
basis.DT<-get_gwas_data(args$manifest_file,args$snp_support_file,args$gwas_data_dir)
#SLE data is a bit ropey and contains p-vals that are zero - not sure why but can
#compute using beta and standard error. I checked and makes little different to PCA's
basis.DT[p.value==0,]$p.value<-1e-13
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.DT[,c('chr','position'):=tstrsplit(pid,':')]
shrink.DT[,c('chr','position'):=tstrsplit(pid,':')]
basis.DT[,position:=as.numeric(position)]

e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

hgnc_name<-'CTLA4'

target<-BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        external_name = hgnc_name,
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T,hgnc_symbol=hgnc_name))

rchr <- gsub("^chr","",unique(seqnames(target)))
rstart<-min(start(range(target)))
rend<-max(end(range(target)))

ideo<-fread("/rds/user/ob219/hpc-work/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

ldb<-unique(subset(basis.DT,chr==rchr & between(position,rstart-1e4,rend+1e4))$ld)
raw.dat<-subset(basis.DT,ld.block %in% c(ldb-1,ldb,ldb+1))
raw.dat[,slp:=sign(log(or)) * -log10(p.value)]
melt.DT<-melt(raw.dat,id.vars=c('pid','chr','position','trait'),measure.vars=c('slp'))
slp<-dcast(melt.DT[variable=="slp",],pid+chr+position~trait+variable)
comb.gr <- with(slp,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(comb.gr)<- as.data.frame(slp[,4:ncol(slp)]) %>% DataFrame
bshrink.gr <- with(shrink.DT[pid %in% slp$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),shrink=bshrink))
wshrink.gr <- with(shrink.DT[pid %in% slp$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),wshrink=ws_ppi))
#compute Wakefields
ppi<-raw.dat[,list(pid=pid,ppi=wakefield_pp(p.value,maf,unique(n),unique(n1/n))),by=c('trait','ld.block')][,.(trait,ppi,pid,ld.block)]
setkeyv(ppi,c('pid','trait'))
setkeyv(raw.dat,c('pid','trait'))
raw.dat <- ppi[raw.dat]
melt.DT<-melt(raw.dat,id.vars=c('pid','chr','position','trait'),measure.vars=c('slp','ppi'))
ppid<-dcast(melt.DT[variable=="ppi",],pid+chr+position~trait+variable)
ppi.gr <- with(ppid,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(ppi.gr)<- as.data.frame(ppid[,4:ncol(ppid)]) %>% DataFrame


tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(comb.gr),bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        chromosome = seqlevels(comb.gr), start = min(slp$position), end = max(slp$position),
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T))
tracks$combmlp<- DataTrack(name="-log10(P)",comb.gr,groups=gsub("_slp","",names(mcols(comb.gr))))
displayPars(tracks$combmlp)$background.title <- 'black'
tracks$ppi <- DataTrack(name="PPA",ppi.gr,groups=gsub("_ppi","",names(mcols(ppi.gr))))
tracks$shrink <- DataTrack(name="PPA Shrink",bshrink.gr,col='dodgerblue')
displayPars(tracks$shrink)$background.title <- 'dodgerblue'
tracks$shrink2 <- DataTrack(name="Weighted PPA Shrink",wshrink.gr,col='firebrick')
displayPars(tracks$shrink2)$background.title <- 'firebrick'
#plotTracks(tracks,from=min(raw.dat[ld.block %in% ldb,]$position),to=max(raw.dat[ld.block %in% ldb,]$position))

## zoom into ld block with highest posterior
tmp.DT<-shrink.DT[ld.block %in% ldb,list(tot=sum(ws_ppi)),by=ld.block]
max.ldb <- tmp.DT[which.max(tot),]$ld.block
pdf("~/tmp/ctla4_shrinkage.pdf")
plotTracks(tracks,from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
dev.off()
