#It is critical for our application that OR are wrt to the same allele
#here is a method for checking this.
library(data.table)
setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")


## first get all of the information that we can from immunobase
traits=c('ATD','CEL','CRO','JIA','MS','PBC','PSC','PSO','RA','SLE','T1D','UC')
curated<-rbindlist(lapply(traits,function(t){
  url<-sprintf("https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-%s-assoc_variantsTAB",t)
  tmp<-fread(url,skip=5L)
  tmp$trait<-t
  md<-names(sort(table(tmp$Disease),decreasing=TRUE)[1])
  message(md)
  subset(tmp,Disease==md)
}))


curated[,c('maj','min'):=tstrsplit(curated$Alleles,'>')]
## the OR are wrt to minor allele so we can compute the risk allele
## remove the instances where the allele is not defined
curated<-curated[!is.na(curated$maj),]
curated[,risk.allele:=min]
setnames(curated,make.names(names(curated)))
curated[curated$Odds.Ratio<1,risk.allele:=maj]
curated[,chr:=sub("([^pq]+)[pq].*","\\1",Region)]
curated[,id:=paste(chr,Position,sep=':')]
curated<-curated[curated$P.Value<5e-8,]


DT<-getGWASData()

DT[DT$disease=='meta.t1d',disease:='T1D']
DT[DT$disease=='CD',disease:='CRO']

DT<-subset(DT,disease %in% unique(curated$trait))

curated<-split(curated,curated$trait)
DTs<-split(DT,DT$disease)

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

disagree<-lapply(seq_along(DTs),function(i){
  dname<-names(DTs)[i]
  message(dname)
  d<-DTs[[i]]
  setkey(d,id)
  cur<-curated[[dname]][,.(id,risk.allele,trait,Odds.Ratio)]
  setkey(cur,id)
  tmp<-subset(d[cur],!is.na(name))
  ## try fixing risk alleles assuming we got it right
  flip.idx<-which(tmp$or<1)
  if(length(flip.idx)>0){
    tmpa<-tmp[flip.idx,]$a1
    tmp[flip.idx,]$a1<-tmp[flip.idx,]$a2
    tmp[flip.idx,]$a2<-tmpa
  }
  tmp<-tmp[tmp$a1 != tmp$i.risk.allele,.(trait,name,a1,a2,risk.allele,risk.allele.freq,maf,i.risk.allele,Odds.Ratio)]
  tmp<-tmp[tmp$a1 != comp(tmp$i.risk.allele),]
})




#Hs_GRCh37-CEL-assoc_variantsTAB
