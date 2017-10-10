library(data.table)
library(GenomicRanges)


#' Compute minor allele frequency shrinkage
#' \code{maf_shrinkage} computes a shrinkage metric for a given list of minor allele frequencies'
#'
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @return a vector of shrinkage metrics
#' @export

maf_shrinkage <- function(f){
  1/sqrt(f * (1-f))
}

#' helper function to sum logs without loss of precision
#' \code{logsum} sums logs without loss of precision
#'
#' @param x a vector of logs to sum
#' @return a scalar

logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

#' compute posterior probabilities using Wakefield's approximate Bayes Factors
#' \code{wakefield_pp} computes posterior probabilities for a given SNP to be causal for a given SNP under the assumption of a single causal variant.
#'
#' @param p a vector of univariate pvalues from a GWAS
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @param N a scalar or vector for total sample size of GWAS
#' @param s a scalar representing the proportion of cases (n.cases/N)
#' @param pi a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
#' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
#' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
#' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
#' is in the range of 0.66-1.5 at any causal variant.
#' @return a vector of posterior probabilities.
#' @export

wakefield_pp <- function(p,f, N, s,pi=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}

#' compute reciprical posterior probabilities using Wakefield's approximate Bayes Factors
#' \code{wakefield_null_pp} computes posterior probabilities for a given SNP to be NOT be causal for a given SNP under the assumption of a single causal variant.
#'
#' @param p a vector of univariate pvalues from a GWAS
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @param N a scalar or vector for total sample size of GWAS
#' @param s a scalar representing the proportion of cases (n.cases/N)
#' @param pi a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
#' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
#' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
#' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
#' is in the range of 0.66-1.5 at any causal variant.
#' @return a vector of posterior probabilities.
#' @export

wakefield_null_pp <- function(p,f, N, s,pi=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    po = exp(lABF + log(p) - log(1-p))
    pp = po/(1+po)
}

#' This function computes the posterior prob that a SNP is causal in a set of traits
#' \code{basis_shrinkage} computes the posterior probability that a SNP is causal across a set of traits
#'
#' @param bf a vector of approximate Bayes Factors using Wakefield's method.
#' @param a scalar or vector of posterior probabilites
#' @return a vector of combined posterior probabilities

basis_shrinkage<-function(BF,pi_i){
  lABF<-log(BF)
  tABF <- c(lABF,0)
  vpi_i<-c(rep(pi_i,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(pi_i)-sBF)
}

#' convert Z score to a signed p value

p2z <- function(p,lor){
  z <- qnorm(0.5 * p.val, lower.tail = FALSE)
  if(missing(lor))
    return(z)
  return(z * sign(lor))
}

#' convert z to p value

z2p <- function(z){
  2* pnorm(abs(z), lower.tail = FALSE)
}

# This function reads a GWAS into a data table

read_GWAS <- function(f,label){
  #if label is missing just use the filename
  if(!file.exist(f))
    stop(sprintf("Cannot find file %s",f))
  DT<-fread(f)
  if(missing(label))
    label<-gsub(".csv","",basename(f),fixed=TRUE)
  DT[,trait:=label]
  return(DT)
}

# this function reads in maf data

add_ref_maf <- function(snp_support_file,DT){
  if(!file.exists(snp_support_file))
    stop(sprintf("Cannot find file %s",snp_support_file))
  ss<-fread(snp_support_file)
  ## use data table to merge the two files
  ss[,pid:=paste(chr,position,sep=':')]
  ss<-ss[,.(pid,maf)]
  setkey(ss,pid)
  tmp<-DT[ss]
  if(nrow(tmp)!=nrow(DT))
    stop("Something went wrong perhaps there are duplicates (by position) in your snp support file or in GWAS input")
  #return(tmp)
}

add_ld_block <- function(ld_support_file,DT){
  if(!file.exists(ld_support_file))
    stop(sprintf("Cannot find file %s",ld_support_file))
  ld<-fread(ld_support_file)
  ld.gr<-with(ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
  snps.gr<-with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,width=1)))
  ol<-as.matrix(findOverlaps(snps.gr,ld.gr))
  DT[ol[,1],ld.block := ol[,2]]
}

# this function gets GWAS data using a manifest file. If a trait list is supplied
# then gets just those traits, if trait list is missing assumes that you want just
# basis

get_gwas_data <- function(manifest_file,snp_manifest_file,ld_support_file,data_dir,trait_list){
  if(missing(trait_list)){
    man<-fread(manifest_file)[basis_trait==1,]
  }else{
    man<-fread(manifest_file)[trait %in% trait_list,]
  }
  man[,file:=file.path(data.dir,file)]
  ret<-rbindlist(lapply(1:nrow(man),function(i){
    message(sprintf("Processing %s",man[i,]$trait))
    tDT<-fread(man[i,]$file)
    tDT[,pid:=paste(chr,position,sep=':')]
    tDT[,c('trait','n','n1') := man[i,.(trait,cases+controls,cases)]]
  }))
  setkey(ret,pid)
  ## next add minor allele frequencies
  message("Adding reference minor allele freq.")
  add_ref_maf(snp_manifest_file,ret)
  message("Assigning LD Blocks")
  add_ld_block(ld_support_file,ret)
  ret
}

## tomorrow use the above routine in a script and then compute the necessary other stats
