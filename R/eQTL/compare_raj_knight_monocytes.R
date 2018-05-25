## let's compare Raj and knight monocyte analysis

## load in knight dataset

library(data.table)
library(magrittr)

raj <- readRDS("/home/ob219/scratch/as_basis/raj_cd14_ind_analysis/regression_analysis/cd14.coeff.pvals.RDS") %>% rbindlist
knight <- readRDS("/home/ob219/scratch/as_basis/knight_cd14_ind_analysis/regression_analysis/cd14.coeff.pvals.RDS") %>% rbindlist


## next get ensg id for each study


kanno <- get(load('/scratch/wallace/twas/knight-probes.RData')) %>% data.table
kanno <- kanno[!is.na(ensembl_gene_id),.(Name,ensembl_gene_id,hgnc_symbol)]
setnames(kanno,c('probe.id','ensg','gene'))
ranno <- readRDS("/scratch/wallace/twas/gene_details.rds") %>% data.table
ranno <-  ranno[,.(affy_hugene_1_0_st_v1,ensembl_gene_id,external_gene_name)]
setnames(ranno,c('probe.id','ensg','gene'))
ranno <- ranno[!duplicated(probe.id),]
kanno <- kanno[!duplicated(probe.id)]
setkey(ranno,probe.id)
setkey(kanno,probe.id)
raj[,probe.id:=as.numeric(gsub("P\\_","",probe))]
setkey(raj,probe.id)
mraj <- ranno[raj][!is.na(ensg),]
setkey(knight,probe)
mknight <- kanno[knight][!is.na(ensg),]

mknight[,id:=paste(PC,ensg,sep=':')]
setkey(mknight,id)
mraj[,id:=paste(PC,ensg,sep=':')]
setkey(mraj,id)
setnames(mraj,'p.coeff','raj.p.coeff')
setnames(mknight,'p.coeff','knight.p.coeff')
o.m<-mknight[mraj]

library(ggplot2)

ggplot(o.m,aes(x=-log10(knight.p.coeff),y=-log10(raj.p.coeff))) + geom_point() + facet_wrap(~PC)

sapply(split(o.m,o.m$PC),function(x) cor(-log10(x$knight.p.coeff),-log10(x$raj.p.coeff)))


pls <- readRDS("/scratch/wallace/twas/probe_alignment.rds") %>% data.table
ranno <- pls[,.(probe,ensembl_gene_id.x,hgnc_symbol)]
setnames(ranno,c('probe.id','ensg','gene'))
kanno <- pls[,.(Name,ensembl_gene_id.x,hgnc_symbol)]
setnames(kanno,c('probe.id','ensg','gene'))
ranno <- ranno[!duplicated(probe.id),]
kanno <- kanno[!duplicated(probe.id)]
setkey(ranno,probe.id)
setkey(kanno,probe.id)


mraj <- ranno[raj][!is.na(ensg),]
mknight <- kanno[knight][!is.na(ensg),]

mknight[,id:=paste(PC,ensg,sep=':')]
setkey(mknight,id)
mraj[,id:=paste(PC,ensg,sep=':')]
setkey(mraj,id)
setnames(mraj,'p.coeff','raj.p.coeff')
setnames(mknight,'p.coeff','knight.p.coeff')
o.m<-mknight[mraj]

library(ggplot2)

ggplot(o.m,aes(x=-log10(knight.p.coeff),y=-log10(raj.p.coeff))) + geom_point() + facet_wrap(~PC)

sapply(split(o.m,o.mOK$PC),function(x) cor(-log10(x$knight.p.coeff),-log10(x$raj.p.coeff)))
