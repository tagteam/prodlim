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

# in a competing risk model multiple types of failure occur
# --------------------------------------------------------------------

compRiskHistory1 <- Hist(time=1:5,event=c("cause 1","cause 2","cause 3", 0,"cause 4"))
plot(compRiskHistory1,layout=list(nrow=3,ncol=3,box.pos=list(c(2,2),c(1,1),c(1,3),c(3,1),c(3,3))),arrow.lab.offset=5)

compRiskHistory2 <- Hist(time=1:5,event=c("Cancer","Heart attack","Car\naccident", 0,"Murder"))
plot(compRiskHistory2,layout=list(nrow=3,ncol=3,box.pos=list(c(2,2),c(1,1),c(1,3),c(3,1),c(3,3))),arrow.lab.offset=5,state.cex=1.8,xbox.rule=.1,ybox.rule=.4,tagBoxes=TRUE)

# in a multi state model events occur in a certain order
# --------------------------------------------------------------------

## illness-death model without recovery
illness.death.frame <- data.frame(time=1:4,
		    from=c("Disease-free","Disease-free",
		      "Diseased","Disease-free"),
		    to=c("0","Diseased","Dead","Dead"))

IDHist <- with(illness.death.frame,Hist(time,event=list(from,to)))
plot(IDHist,ybox.rule=4,xbox.rule=.3,state.cex=1.3,enum=TRUE,arrow.lab.side=c(-1,-1,1))


## illness-death with recovery
illness.death.frame2 <- data.frame(time=1:5,
				   from=c("Disease\nfree","Disease\nfree",
                                     "Diseased","Diseased","Disease\nfree"),
				   to=c("0","Diseased","Disease\nfree",
                                                           "Dead","Dead"))
IDHist2 <- with(illness.death.frame2,Hist(time,event=list(from,to)))
plot(IDHist2,
     ybox.rule=1.3,
     xbox.rule=.3,
     state.cex=2,
     arrow.lab.offset=c(13,13,8,10),
     enum=TRUE,
     verbose=FALSE)

## change the layout of the graphic

plot(IDHist2,
     ybox.rule=1.3,
     xbox.rule=.3,
     state.cex=2,
     enum=TRUE,
     verbose=FALSE,
     layout=list(ncol=3,nrow=2,box.pos=list(c(1,1),c(2,2),c(1,3))),
     arrow.lab.side=c(-1,1,-1,1),
     arrow.lab.offset=c(15,15,10,10))


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

