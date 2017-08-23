library(ggplot2)
devtools::install_github('stefano-meschiari/latex2exp')
library(latex2exp)
## here we examine the distribution of logit and normal quantiles.
x<-c(exp(seq(-30,0.01,1)),seq(from=0.01,to=1,by=0.01))
logit<-function(p) log(p/(1-p))
df<-data.frame(y=logit(x),x=qnorm(x))
df<-df[!is.infinite(df$x),]

ggplot(data = df, aes(x = x, y = y)) +
    geom_line() + 
    geom_abline(slope=1.81,color='red') + theme_bw() + xlab("Normal Quantiles") + ylab("Logit") +
    geom_text(x = -5, y = -5, label = TeX('$\\frac{\\pi}{\\sqrt{3}}$',output="character"),color='red',parse = TRUE) + geom_point(y=logit(5e-8),x=qnorm(5e-8))

## what about the variance of 
