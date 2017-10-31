## find a, b parameters for prior of cases, given E(lor)=0, P(abs(lor)>log(target.or)) = target.prob
est.a1b1 <- function(target.or, target.prob, a, b) {
    ## E.lor given a1,b1
    fn.b <- function(b1,a1) {
        abs(digamma(a) - digamma(b) - digamma(a1) + digamma(b1))
    }

    p0 <- rbeta(10000,shape1=a,shape2=b)
    fn <- function(a1) {
        ## calc b1 for given a1
        b1 <- optimise(fn.b,c(1,b),a1=a1)$minimum
        p1 <- rbeta(10000,shape1=a1,shape2=b1)
        lor <- log(p0) - log(1-p0) + log(1-p1) - log(p1)
        ret <- mean(abs(lor) > log(target.or))
        ## points(a1,ret)
        ret
    }
    f <- a/(a + b)
    a1 <- 1
    ## plot(c(1,1000),c(0,1),type="n"); box(); axis(1:2); abline(h=target.prob)
    while(fn(a1) > target.prob)
        a1 <- a1+1
    c(a1=a1, b1=optimise(fn.b,c(1,b),a1=a1)$minimum)
}


lor.f <- function(f,nchr=5000,nsim=10000) {
    v=f*(1-f)/nchr

    ## prior in controls
    a=-f*(f^2-f+v)/v
    b=(f^2-f+v)*(f-1)/v
    p0 <- rbeta(nsim,shape1=a,shape2=b) ## check - should have mean around f

    ## prior in cases
    a1b1 <- est.a1b1(target.or=2,target.prob = 0.01,a=a,b=b)
    a1 <- a1b1[["a1"]]
    b1 <- a1b1[["b1"]]
    p1 <- rbeta(nsim,shape1=a1,shape2=b1)

    ## posterior in cases if see 00
    sim.lor <- function(n2) {
        posta=2-n2 + a1
        postb=n2 + b1
        p1 <- rbeta(nsim,shape1=posta,shape2=postb)
        ## p1 <- rbeta(nsim,shape1=a1,shape2=b1)
        log(p0*(1-p1)/(p1*(1-p0)))
    }
    lor.00 <- sim.lor(0)
    lor.01 <- sim.lor(1)
    lor.11 <- sim.lor(2)

    ## therefore, posterior expected log or will be
    c("00"=mean(lor.00),
      "01"=mean(lor.01),
      "11"=mean(lor.11),
      ## quantile(1-p1,c(0.05,0.25,0.5,0.75,0.95)),
      quantile(p1,c(0.05,0.25,0.5,0.75,0.95)))
}

af <- seq(0.01,0.99,by=0.01)
library(parallel)
options(mc.cores=6)
results <- as.data.frame(cbind(f=af,do.call("rbind",mclapply(af, lor.f, nsim=100000))))
colnames(results) <- make.unique(colnames(results))
library(ggplot2)
m <- reshape2::melt(results,"f")
head(m)

ggplot(m[grep("%", m$variable,invert=TRUE),],
       aes(x=f,y=value,col=variable,group=variable)) + geom_point() + geom_smooth() + labs(x="AF",y="log OR")

ggplot(m[grep("%", m$variable),],
       aes(x=f,y=value,col=variable,group=variable)) + geom_path() + labs(x="ctl AF",y="case AF") +
  geom_abline()
