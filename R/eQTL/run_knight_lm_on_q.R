library(data.table)


## this code is too slow and get's me banned if I run interactively so let's submit to
##q

for.reg <- readRDS("/home/ob219/scratch/as_basis/knight_cd14_ind_analysis/regression_analysis/cd14.RDS")


## first we attempt to fit a set of linear models where PC is the explantory model


## for the number of probes this code becomes unwieldy and cannot be run interactively

all.lms <- lapply(names(for.reg)[grep("PC[0-9]+",names(for.reg))],function(pc){
  message(pc)
  lapply(grep('^ILMN\\_[0-9]+',names(for.reg)),function(i){
    probe <- names(for.reg)[i]
    dat <- for.reg[,c(pc,probe),with=FALSE]
    setnames(dat,c('PC','expression'))
    lm(expression~PC,dat)
  })
})

names(all.lms)<-names(for.reg)[grep("PC",names(for.reg))]
## for each model we can obtain t-statistics for the beta coefficients as p-values.
## we can the check for inflation by plotting as a qqplot trellis plot.

## should create a library for this.
get_qq_dt <- function(x,l=0.99,minx=TRUE){
  n <- length(x)
  ## expected
  q <- -log10((n:1)/(n+1))
  ## observed
  x <- sort(x)
  ## l is a parameter that allows us to only plot a subset of the data
  n1 <- round(l*n)
  if(minx)
    return(data.table(expected=c(0,q[n1:n]),observed=c(0,x[n1:n])))
  return(data.table(expected=q,observed=x))
}

pc.p <- lapply(seq_along(names(all.lms)),function(i){
  message(i)
  pc<-names(all.lms)[i]
  all.p <- sapply(all.lms[[i]],function(x){
    summary(x)$coefficient["PC",4]
  })
  DT <- data.table(PC=pc,probe=names(for.reg)[grep('^ILMN\\_',names(for.reg))],p.coeff=all.p)
})

saveRDS(pc.p,file="/home/ob219/scratch/as_basis/knight_cd14_ind_analysis/regression_analysis/cd14.coeff.pvals.RDS")
