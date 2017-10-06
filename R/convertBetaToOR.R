## code to convert linear regression summary stats to their logistic counterparts

# compute counts for different alleles within a population
# N - Sample size
# m - Minor allele frequency
n.hom<-function(N,m) N * m^2
n.het<-function(N,m) 2 * N * m * (1-m)
n.wt<-function(N,m) N * (1-m)^2

## estimate Y mean
## b - estimate of beta
# N - Sample size
# m - Minor allele frequency

y_bar<-function(b,N,m) (b * (n.het(N,m) + (2 * n.hom(N,m))))/N

## estimate allelic sum of squares
# N - Sample size
# m - Minor allele frequency

ssx<-function(N,m) 2 * m * (N-1) * (1-m)

# standard error of epsilon
# N - Sample size
# m - Minor allele frequency
# seb - Standard error for beta
sigma_e<-function(N,seb,m) seb * sqrt(ssx(N,m))

## estimate P(Y>thresh|X=x)
# N - Sample size
# m - Minor allele frequency
# seb - Standard error for beta
# thresh - threshold to use
threshP<-function(thresh,X,N,b,seb,m){
  SE.epsilon <- sigma_e(N,seb,m)
  pnorm(thresh,mean=X*b,sd=SE.epsilon,lower.tail = FALSE)
}

## estimate P(Y>thresh|X=x)
# N - Sample size
# m - Minor allele frequency
# seb - Standard error for beta
# thresh
threshOdds<-function(thresh,X,N,b,seb,m){
  num<-threshP(thresh,X,N,b,seb,m)
  num/(1-num)
}

## compute log(OR) for threshold dichotomised trait
# N - Sample size
# m - Minor allele frequency
# seb - Standard error for beta
# thresh - threshold to use if missing uses ybar
threshBeta<-function(thresh,N,b,seb,m){
  lyx1<-log(threshOdds(thresh,1,N,b,seb,m))
  lyx0<-log(threshOdds(thresh,0,N,b,seb,m))
  lyx1-lyx0
}


## compute the stderr of thresh beta
threshSigmaBeta<-function(thresh,n1,n0,m,tb){
  ## number above thresh
  N<-n0+n1
  ploidy<-2 # given current framework this is constant
  a<-(n0 * (1-m))/N
  b<-(n0 * m)/N
  c<-((n1*a)/(a+(b*exp(tb))))/N
  d<-((n1*b)/(a+(b*exp(tb))))/N
  var_maf<-sqrt(rowSums(do.call('cbind',lapply(list(a,b,c,d),function(fi) 1/fi))))
  var_ss<-sqrt(1/N)
  var_ploidy<-sqrt(1/ploidy)
  var_maf * var_ss * var_ploidy
}

## overall conversion routine
convertBetaToOR<-function(thresh,N,b,seb,m){
  if (missing(thresh))
    thresh <- y_bar(b,N,m)
  tb<-threshBeta(thresh,N,b,seb,m)
  x_0<-threshP(thresh,0,N,b,seb,m) * n.wt(N,m)
  x_1<-threshP(thresh,1,N,b,seb,m) * n.het(N,m)
  x_2<-threshP(thresh,2,N,b,seb,m) * n.hom(N,m)
  n1<-round(x_0 + x_1 + x_2)
  n0<-N-n1
  tsb<-threshSigmaBeta(thresh,n1,n0,m,tb)
  Z<-tb/tsb
  P<-2*(pnorm(abs(Z),lower.tail = FALSE))
  list(OR=exp(tb),beta=tb,se.beta=tsb,P=P,Z=Z,thresh=thresh,n0=n0,n1=n1)
}

## test to see if it works
if(FALSE){
  sN<-10000
  smaf<-0.1
  sbeta<-0.4
  sp.val<-5e-8
  sbeta.se<-sbeta/abs(qnorm(sp.val/2))
  convertBetaToOR(N=sN,b=sbeta,seb=sbeta.se,m=smaf)
}
