library(data.table)
library(ggplot2)
library(ggplot2)


SIM.DIR<-'/home/ob219/scratch/as_basis/tmp/sample_size/'

fs<-list.files(path=SIM.DIR,pattern='*.RData',full.names=TRUE)

res<-lapply(fs,function(f){
  dat<-get(load(f))
  rbind(data.table(metric='gamma_hat_sample_size',case.no=dat$case.no,control.no=dat$ctrl.no,simulation=dat$sims$unscaled_gh_ss),
  data.table(metric='gamma_hat_maf',case.no=dat$case.no,control.no=dat$ctrl.no,simulation=dat$sims$unscaled_gh_maf))
})

res<-rbindlist(res)

pdf("/home/ob219/git/as_basis/pdf/vary_case.pdf")
ggplot(res[res$control.no!=4000,],aes(x=as.factor(case.no),y=simulation)) + geom_boxplot() + facet_wrap(~metric) + theme_bw() + xlab("Number of Cases (Controls = 3000)") + ylab("Simulated ill.t1d eucledian distance")
dev.off()

## can use an empriccal Bayes technique to work out what our prior should be ?
library(reshape2)
library(data.table)
library(matrixStats)


pi_1<-1e-4
setwd("/home/ob219/git/as_basis/R/")
source("./central_functions.R")
source("./wakefield.R")



DT<-getGWASData()
## remove quantitative traits as well as meta.t1d
QT<-c('eosinophil','lymphocyte','myeloid','wbc')
DT<-DT[!DT$disease %in% c(QT,'ill.t1d','aff.t1d'),]
## compute the different Metrics proposed
## beta == log(or) DT$lor already present
## Just use Z score
DT[,Z:=qnorm(0.5 * p.val, lower.tail = FALSE) * sign(lor)]
DT[,lp0:=log(1-approx.bf.z2(Z,maf,cases+controls,cases/(cases+controls),pi_1)),by=c('disease','ld.block')]

## next create a matrix that we can sample from

tmp<-melt(DT,id.vars=c('id','disease'),measure.vars = 'lp0')
mat<-dcast(tmp,id~disease+variable)
mat$id<-NULL
## OK next sample different collections of diseases and compute q_i

mat<-as.matrix(mat)

eb<-lapply(1:length(colnames(mat)),function(nd){
     message(sprintf("Checking %d",nd))
    if(nd==1)
      return(1-exp(mat))
    ## get the combination matrix
    cmat<-combn(1:length(colnames(mat)),nd)
    message(sprintf("Checking %d with %d combinations",nd,ncol(cmat)))
    apply(cmat,2,function(i) {
      1-exp(rowSums(mat[,i]))
    })
})

for.plot<-do.call('rbind',lapply(seq_along(eb),function(i){
  cbind(i,colMeans(eb[[i]]))
}))

for.plot<-as.data.table(for.plot)
setnames(for.plot,c('n','q'))

ggplot(for.plot,aes(x=n,y=q,z=n2,p=n3)) + geom_point() + xlab("Number of diseases") + ylab(expression(bar(q[i])))


## code for computing the best fit for prior computation

for.plot$n2<-for.plot$n^2
for.plot$n3<-for.plot$n^3

## here we work out the best model for predicting pi based on number of diseases

mod1<-lm(q~n-1,data=for.plot)
mod2<-lm(q~n + n2 -1,data=for.plot)
mod3<-lm(q~n + n2 + n3 -1,data=for.plot)

lapply(list(mod1,mod2,mod3),summary)

## Likelihood of data given model always improves as you add more terms. BIC (like AIC) penalise the increase by an amount that depends on the size of the data. Smallest BIC wins. Allows comparison of non nested models, unlike ANOVA
lapply(list(mod1,mod2,mod3),BIC)

## this means that model two wins i.e q = (6.6e-4 x n) - (1.3 x n^2)
