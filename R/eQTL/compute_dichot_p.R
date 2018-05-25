library(data.table)
library(magrittr)

## this is some code to compute empirical p.values for dichot trait PC loadings
## convention is to put actual data in 'projections' and NULL data in 'projection_nullX'

# dir containing the projections
pdir <- '/home/ob219/scratch/as_basis/raj_cd4_dichot_analysis/'
# dir containing eigen basis for AI (learned from same SNPs as used to project dichotomous trait)
pc.file=file.path(pdir,"pc.RDS")
pc.emp <- readRDS(pc.file)
act.dir <- file.path(pdir,"projections")
act.proj <- list.files(path=act.dir,pattern="*.RDS",full.names=TRUE) %>% lapply(.,readRDS) %>%
          do.call('rbind',.)
## next get null
null.dirs <- list.dirs(path = pdir, full.names = TRUE) %>% .[grep("projections_null",.)]
## read these in
null.proj <- sapply(null.dirs,list.files,full.names=TRUE,pattern="*.RDS") %>% lapply(.,readRDS) %>%
          do.call('rbind',.)
## next we want to generate a p.value for each probe and principal component loading


## get vector of control
ctrl<-pc.emp$x["control",]

compD <- function(M,cD){
  M - matrix(rep(cD,times=nrow(M)),ncol=length(cD),byrow=TRUE)

}

null.d <- compD(null.proj,ctrl)
act.d <- compD(act.proj,ctrl)

p.val <- act.d
p.val <- 0

for(i in 1:ncol(act.d)){
  message(i)
  p.val[,i] <- sum(act.d[,i] > null.d[,i]) / nrow(null.d)
  idx <- which(act.d[,i]<0)
  p.val[idx,i] <- 1 - p.val[idx,i]
}
