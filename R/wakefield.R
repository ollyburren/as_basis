logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

logdiff <- function(x,y) {
    my.max <- max(x,y)
    my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
    return(my.res)
}

## compute variance shrinkage for quantitative trait study
Var.data <- function(f, N) {
    1 / (2 * N * f * (1 - f))
}

## compute variance shrinkage for case control study
Var.data.cc <- function(f, N, s) {
    1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(p,f, N, s,pi_i,type='CC') {
    if(type=="QUANT") {
      sd.prior <- 0.15
      V <- Var.data(f, N)
    } else {
      sd.prior <- 0.2
      V <- Var.data.cc(f, N, s)
    }
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    sBF <- logsum(lABF + log(pi_i))
    exp(lABF + log(pi_i))/(exp(sBF) + 1)
}

approx.bf.z <- function(z,f, N, s,pi_i) {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}
