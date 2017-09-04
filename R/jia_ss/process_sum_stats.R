library(data.table)

in.dir<-'/scratch/ob219/jia/summary-stats/'
out.dir<-'/scratch/ob219/jia/by.trait/'

fs<-list.files(path=in.dir,pattern='*.RDS',full.names=TRUE)
by.an<-split(fs,gsub(".*[_]([^\\.]+).*","\\1",basename(fs)))

processTrait<-function(tl){
  tmp<-rbindlist(lapply(tl, function(f){
  message(sprintf("Reading %s",f))
  readRDS(f)
  }))
  tmp[order(tmp$chr,tmp$position),]
}


for(i in seq_along(by.an)){
  DT<-processTrait(by.an[[i]])
  DT<-DT[,.(rsid,chr,position,a1,a2,alt,af.wrt.alt,beta,Var.beta,Z,p.val)]
  saveRDS(DT,file=file.path(out.dir,paste(names(by.an)[i],'RDS',sep='.')))
}
