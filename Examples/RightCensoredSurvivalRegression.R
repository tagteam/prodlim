library(SimSurv)
library(prodlim)
set.seed(8)
SurvFrame <- SimSurv(N=100,
                     surv=list(dist = "rweibull",
                       args = list(shape=1),
                       baseline = 1/100,
                       link="exp",
                       coef=1),
                     cens=list(dist = "rexp",
                       args=NULL,
                       baseline=1/1000,
                       link="exp",
                       coef=0,
                       max = 150,
                       type = "interval",
                       lateness=1,
                       unit=10),
                     cova=list(X1 = list("rnorm", mean = 0, sd = 2),
                       X2 = list("rbinom", size=1,prob=.5)),
                     keep.uncensored=TRUE,
                     method=c("simulation"),
                     verbose=1)
 print(SurvFrame,digits=2)

fitPL <- prodlim(formula=Hist(time=list(L,R),event=status!=0)~1,
        data=SurvFrame,
        maxiter=1,
##         grid=,
        tol=7,
        ml=FALSE)

fitTRUE <- prodlim(Hist(time=uncensored.time,
                        event=rep(1,NROW(SurvFrame)))~1,data=SurvFrame)

fitML <- prodlim(formula=Hist(time=list(L,R),event=status!=0)~1,
                 data=SurvFrame,
                 maxiter=1000,
                 ##         grid=,
                 tol=7,
                 ml=TRUE)

plot(fitML,col=3)
plot(fitPL,add=T,col=2)
plot(fitTRUE,add=T,col=1)

