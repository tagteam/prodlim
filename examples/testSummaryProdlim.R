library(prodlim)
library(survival)

# {{{survival data

# {{{ no covariates

d <- data.frame(time=c(1,2,2,3,4,8,8,8,9),status=c(0,1,0,1,0,1,0,1,1))
f <- prodlim(Hist(time,status)~1,data=d)
s <- summary(survfit(Surv(time,status)~1,data=d))
y <- cbind(s$time,s$n.risk,s$n.event,s$n.censor,s$surv,s$std.err,s$lower,s$upper)
x <- summary(f,times=s$time)
all.equal(y,x$table)

x <- summary(f,times=s$time)$table
summary(f,intervals=TRUE,times=c(0,8,9,10))
stopifnot(all(data.frame()))
# }}}
# {{{ with covariates
d <- SimSurv(100)
summary(prodlim(Hist(time,status)~X1,data=d))
summary(prodlim(Hist(time,status)~X1,data=d),newdata=data.frame(X1=c(1)))
# }}}
# {{{ clustered data 
set.seed(17)
d <- SimSurv(10)
d$id <-  rbinom(10,6,.5)
f <- prodlim(Hist(time,status)~X2,data=d)
if (!is.function("cluster")) cluster <- function(x)x
g <- prodlim(Hist(time,status)~cluster(id),data=d)
fg <- prodlim(Hist(time,status)~X2+cluster(id),data=d)
summary(fg)
# }}}
# }}}

# {{{competing risk data

# {{{ no covariates
D <- SimCompRisk(100)
summary(prodlim(Hist(time,event)~1,data=D))
summary(prodlim(Hist(time,event)~1,data=D),cause=1)
summary(prodlim(Hist(time,event)~1,data=D),cause=1:2)
summary(prodlim(Hist(time,event)~1,data=D),cause="2")
# }}}

# {{{ with covariates
D <- SimCompRisk(100)
D$E <- factor(D$event,levels=c(0,1,2),labels=c("0","a","b"))
summary(prodlim(Hist(time,event)~X1,data=D))
summary(prodlim(Hist(time,event)~X1,data=D),cause=1:2,newdata=data.frame(X1=c(0,1)))
summary(prodlim(Hist(time,event)~X1,data=D),cause=1:2,newdata=data.frame(X1=c(1)))
summary(prodlim(Hist(time,event)~X1,data=D),cause=1,newdata=data.frame(X1=c(1:0)))
summary(prodlim(Hist(time,event)~X1,data=D),cause=1,newdata=data.frame(X1=c(1:0)))
summary(prodlim(Hist(time,E)~X1,data=D),cause=1,newdata=data.frame(X1=c(1)))

summary(prodlim(Hist(time,event)~X2,data=D))
summary(prodlim(Hist(time,event)~X2,data=D),cause=1:2,newdata=data.frame(X2=c(0,1)))
summary(prodlim(Hist(time,event)~X2,data=D),cause=1:2,newdata=data.frame(X2=c(1)))
summary(prodlim(Hist(time,event)~X2,data=D),cause=1,newdata=data.frame(X2=c(1:0)))
summary(prodlim(Hist(time,event)~X2,data=D),cause=1,newdata=data.frame(X2=c(0,1,0)))
summary(prodlim(Hist(time,E)~X2,data=D),cause=1,newdata=data.frame(X2=c(1)))

summary(prodlim(Hist(time,event)~X1+X2,data=D))
summary(prodlim(Hist(time,event)~X1+X2,data=D),cause=1:2,newdata=data.frame(X1=c(0,1),X2=0.1))
summary(prodlim(Hist(time,event)~X1+X2,data=D),cause=1:2,newdata=data.frame(X1=c(1),X2=c(0.1,0.4)))
summary(prodlim(Hist(time,event)~X1+X2,data=D),cause=1,newdata=data.frame(X1=c(1),X2=c(0.1,0.4)))
summary(prodlim(Hist(time,event)~X1+X2,data=D),cause=1,newdata=data.frame(X1=c(1),X2=c(0.1,0.4)))
summary(prodlim(Hist(time,E)~X1+X2,data=D),cause=1,newdata=data.frame(X1=c(0,1,0),X2=c(0.1,0.4,0.8)),times=17)
summary(prodlim(Hist(time,E)~X1+X2,data=D),cause=1,newdata=data.frame(X1=c(0,1,0),X2=c(0.1,0.4,0.8)),times=c(17,0,100))
# }}}

# }}}
