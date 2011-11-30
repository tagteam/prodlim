library(prodlim)
library(survival)

# {{{survival data

# {{{ no covariates

d <- data.frame(time=c(1,2,2,3,4,8,8,8,9),status=c(0,1,0,1,0,1,0,1,1))
f <- prodlim(Hist(time,status)~1,data=d)
s <- summary(survfit(Surv(time,status)~1,data=d))
y <- cbind(s$time,s$n.risk,s$n.event,s$n.censor,s$surv,s$std.err,s$lower,s$lower)
x <- summary(f,times=s$time)
all.equal(y,x$table)

x <- summary(f,times=s$time)$table
summary(f,intervals=TRUE,times=c(0,8,9,10))
stopifnot(all(data.frame()))

# }}}

# {{{ clustered data 
set.seed(17)
d <- SimSurv(10)
d$id <-  rbinom(10,6,.5)
f <- prodlim(Hist(time,status)~X2,data=d)
g <- prodlim(Hist(time,status)~cluster(id),data=d)
fg <- prodlim(Hist(time,status)~X2+cluster(id),data=d)
summary(fg)
# }}}
# }}}
