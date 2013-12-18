library(prodlim)

set.seed(1)
a <- sort(SimSurv(10,list(dist="rweibull"),cova=NULL)$uncensored.time)
set.seed(1)
b <- sort(rweibull(10,shape=1))
cbind(a,b)


set.seed(1)
a <- sort(SimSurv(10,surv=list(dist="rweibull",args=list(shape=2)),cova=NULL)$uncensored.time)
set.seed(1)
b <- sort(rweibull(10,shape=2))
cbind(a,b)

d <- SimSurv(100,cova=list(X=list("rbinom",size=1,prob=.5)),surv=list(dist="rweibull",args=list(shape=.2),coef=c(1)))

plot(d)

d1 <- SimSurv(100,surv=list(dist="rweibull",args=list(shape=2)),cova=NULL)
d2 <- SimSurv(100,surv=list(dist="rweibull",args=list(shape=.2)),cova=NULL)
plot(d1)
plot(d2,add=T)
