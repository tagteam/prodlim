## 1.1. survival uncensored
h <- Hist(1:10)
stopifnot(attr(h,"cens.type")=="uncensored")
stopifnot(attr(h,"entry.type")==NULL)
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("time","status"))

## 1.2. survival right censored
h <- Hist(time=1:10,event=c(0,1,0,0,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="uncensored")
stopifnot(attr(h,"entry.type")==NULL)
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("time","status"))

## 1.3. survival right censored and left-truncated
h <- Hist(entry = 0:9, time=1:10,event=c(0,1,0,0,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="right-censored")
stopifnot(attr(h,"entry.type")=="left-truncated")
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("entry", "time","status"))

## 1.4. survival: right censored and left-truncated
h <- Hist(entry = 0:9, time=1:10,event=c(0,1,0,0,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="right-censored")
stopifnot(attr(h,"entry.type")=="left-truncated")
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("entry", "time","status"))

## 1.5. survival: interval censored and right censored and left-truncated
h <- Hist(entry = 0:9, time=list(1:10,b=2:11),event=c(0,1,0,1,0,0,0,0,0,1))
stopifnot(attr(h,"cens.type")=="interval-censored")
stopifnot(attr(h,"entry.type")=="left-truncated")
stopifnot(attr(h,"model")=="survival")
stopifnot(colnames(h)==c("entry", "L","R","status"))


## The function `Hist' (short for history) if a version of the
## the function `Surv' which is defined in the survival package
##
## `Hist' serves `prodlim' in a similar way as `Surv' serves
## `survfit'
##
## `Hist' provides functionality suitable for competing risk
## and other multi state models.
##
## Interval censored observations are not treated in the
## same way as in the survival package.
## 
## the class `Hist' has summary and plot methods

# in a two state survival model 
# --------------------------------------------------------------------

## If there is only one cause of failure, then `Surv' and
## `Hist' do more or less the same. The attributes and the
## print/summary functions are different. 

library(survival)
library(prodlim)
d <- SimSurv(10)
survHistory1 <- Surv(d$time,d$status)
survHistory2 <- Hist(d$time,d$status)

# both objects hold the same information:
stopifnot(all(unclass(survHistory1)==unclass(survHistory2)))

summary(survHistory1)
summary(survHistory2)
plot(survHistory2)

# `Hist' allows any code for right censored

Hist(time=1:2,event=c("Failure","censCode"),cens.code="censCode")

Hist(time=1:2,event=c(1,9),cens.code=9)

# interval censoring
# --------------------------------------------------------------------

# interval censored times can be specified in two ways:

# 1. time is a list with the left and right interval limits
Hist(time=list(1:5,c(1,5,3,4,8)))

# 2. time is a matrix or data.frame with the left and
# right interval limits

Hist(time=matrix(c(1:5,c(1,5,3,4,8)),ncol=2))

Hist(time=data.frame(L=1:5,R=c(1,5,3,4,8)))


# Right censored observations are defined either
# when the right interval limit is missing or infinitive

Hist(time=list(1:5,c(1,5,3,NA,8)))
Hist(time=list(1:5,c(1,Inf,3,4,8)))

# or via `cens.code' and an event indicator

Hist(time=list(1:5,c(1,5,3,4,8)),event=c(1,0,2,0,1),cens.code=0)

#gives the same as

Hist(time=list(1:5,c(1,NA,3,4,8)),event=c(1,0,2,0,1),cens.code=0)
