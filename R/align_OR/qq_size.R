library(data.table)

#Size callibration qq plots

## see how normal our estimates are for blood
all.sims <- list.files(path='/home/ob219/scratch/as_basis/blood_summ_analysis/',pattern='*.RDS',full.name=TRUE)
all.sims <- split(all.sims,sapply(strsplit(basename(all.sims),'_'),"[[",1))

all.res <- lapply(names(all.sims),function(n){
  do.call('rbind',lapply(all.sims[[n]],function(x) data.table(readRDS(x))))
})

library(ggplot2)
library(cowplot)

# for an empirically generated null projection is it normally distributed

res <- all.res[[1]]

library(magrittr)

forp <- lapply(colnames(res),function(x){
  tmp<-res[[x]]
  z <- (tmp-mean(tmp))/sd(tmp)
  t2 <- qqnorm(z,plot.it=FALSE)
  data.table(x=t2$x,y=t2$y,pc=x)
}) %>% rbindlist

p.vals <- lapply(colnames(res),function(x){
  tmp<-res[[x]]
  z <- (tmp-mean(tmp))/sd(tmp)
  pv <- sprintf("P=%s",signif(ks.test(z,"pnorm")$p.value,digits=3))
  data.table(x=-1,y=2,p.value=pv,pc=x)
}) %>% rbindlist

forp[,pc:=factor(pc,levels=unique(forp$pc))]



pm <- ggplot(forp[!pc %in% c('PC10'),],aes(x=x,y=y)) + geom_point(size=0.5,alpha=0.5) +
geom_abline(col='firebrick') + facet_wrap(~pc,nrow=3,ncol=3) + xlab("Expected Z") +
ylab("Observed Z") + geom_text(data = p.vals[!pc %in% c('PC10'),],aes(label = p.value))
save_plot("~/tmp/sample_size_qq.pdf",pm)


## If we simulate how many simulations should we do. Here this indicates that 500 gives us
## a coefficient of variation of 3%
samp.sizes <- c(200,500,1000,2000,10000,1e5)
sres <- lapply(samp.sizes,function(N){
  tmp<-sapply(1:1000,function(i){
    sd(rnorm(N,sd=0.1))
  })
  sd(tmp)/mean(tmp) * 100
})
names(sres) <- samp.sizes
