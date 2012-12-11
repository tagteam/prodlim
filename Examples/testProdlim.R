## The function `Hist' (short for history) generalizes
## the function `Surv' which is defined in the survival package

# first test: a two state survival model 
# --------------------------------------------------------------------

## If the status variable takes only two values, then `Surv' and
## `Hist' do more or less the same. The attributes and the
## print/summary functions are different. 

require(prodlim)
data(pbc)

sH1 <- Hist(pbc$time,pbc$status)
sH2 <- Surv(pbc$time,pbc$status)
stopifnot(all(unclass(sH1)==unclass(sH2)))
summary(sH1)
summary(sH2)

# second test: a competing risk model 
# --------------------------------------------------------------------

require(prodlim)
require(cmprsk)

x <- 1:15
y <- c(2,1,0,1,2,1,0,1,1,0,1,1,2,1,1)
A <- cuminc(x,y)

Hist(x,y)
B <- prodlim(Hist(x,y)~1)
#plot(B)
#plot(A)
#plot(B,cause=c(1,2))

all(A$"1 1"$est == predict(B,times=c(0,A$"1 1"$time[-length(A$"1 1"$time)]))[[1]])
all(A$"1 2"$est == predict(B,cause=2,times=c(0,A$"1 2"$time[-length(A$"1 2"$time)]))[[1]])

test <- data.frame(time=rep(1:10,2),status=rep(0:1,10))
pfit <- prodlim(Hist(time,status)~1,data=test)
summary(pfit,times=c(2,4.4,9))


prodlim(Surv(time,status)~1,data=test) ## same as pfit
fit <- survfit(Surv(time,status)~1,data=test) ## same as pfit$surv
sfit <- summary(fit,times=c(2,4.4,9))
all(round(pfit$surv,10)==round(survfit(Surv(time,status)~1,data=test)$surv,10))

# survival
## require(survival)
## data(pbc)
pfit <- prodlim(Hist(time,status)~1,data=pbc)
prodlim(Surv(time,status)~1,data=pbc) ## same as pfit
survfit(Surv(time,status)~1,data=pbc)$surv ## same as pfit$surv
all(round(pfit$surv,10)==round(survfit(Surv(time,status)~1,data=pbc)$surv,10))
## print(pfit)
## summary(pfit)
## plot(pfit,atrisk=T,conf.int=T)

a <- prodlim(Surv(time,status),data=pbc)
a <- prodlim(Hist(time,status)~1,data=pbc)
b <- prodlim(Surv(time,status)~1,data=pbc,reverse=TRUE)
c <- prodlim(Surv(time,status)~edema,data=pbc,reverse=TRUE)
d <- prodlim(Surv(time,status)~NN(age),data=pbc,reverse=TRUE)

## class(c) <- "dynamic"


normal.cox <- coxph(Surv(time,status)~edema,data=pbc)

update.cox <- function(object,tstar,data){
  object$call$data <- data[data$time>tstar,]
  update <- eval(object$call)
  class(update) <- "dynamicCox"
  update
}

predictProb.dynamicCox <- function(object,newdata,cutpoints,learn.data,...){
  p <- matrix(1,nrow=NROW(newdata),ncol=length(cutpoints))
  p
}


pec.c <- pec(object=list(c),
             formula=Surv(time,status)~1,
             data=pbc,
             exact=TRUE,
             method="ipcw",
             cens.model="marg",
             B=0,
             verbose=TRUE)


a
b
c
d

par(mfrow=c(2,2))
plot(a)
plot(b)
plot(c)
plot(d)


#cluster
if (!is.function("cluster")) cluster <- function(x)x
testdat <- data.frame(time=c(1,2,2,3,7),status=c(0,1,1,0,1),patnr=c(1,2,1,3,2))
y <- prodlim(Hist(time,status)~cluster(patnr),data=testdat)

library(dental)
x <- clustersurv(Surv(time,status)~1|patnr,data=testdat)
x <- prodlim(Hist(survtime,survstatus)~cluster(patnr),data=Sterioss)
data(primary)
data(Sterioss)


w <- as.numeric(Sterioss$patnr)


# bootstrapping study units
# --------------------------------------------------------------------
N <- NROW(Sterioss)
boot <- matrix(sapply(1:B,function(b){sample(1:N,replace=TRUE)}),nrow=N,ncol=B)
bootx <- sapply(1:B,function(b){
  if ((b/100==round(b/100)))print(b)
  u <- predict(prodlim(Hist(survtime,survstatus)~1,data=Sterioss[boot[,b],]),times=seq(0,72,12),verbose=F)$surv
})
# bootstrapping patients
# --------------------------------------------------------------------
B <- 1000
N2 <- length(unique(Sterioss$patnr))
Sterioss.patlist <- split(Sterioss[,c("survstatus","survtime","patnr")],Sterioss$patnr)
boot2 <- matrix(sapply(1:B,function(b){sample(1:N2,replace=TRUE)}),nrow=N2,ncol=B)
bootx2 <- sapply(1:B,function(b){
  if ((b/100==round(b/100)))print(b)
  data2 <- do.call("rbind",Sterioss.patlist[boot2[,b]])
  if(sum(data2$survstatus)==0) print(paste(b,"failed"))
  else predict(prodlim(Hist(survtime,survstatus)~1,data=data2),times=seq(0,72,12),verbose=F)$surv
})

fit <- prodlim(Surv(survtime,survstatus)~cluster(patnr),data=Sterioss)

test <- bootCluster.product.limit(fit,
                                  times=sort(unique(Sterioss$survtime)),
                                  B=10)


# repeated sampling one unit per patient 
# --------------------------------------------------------------------
B <- 1000
bootx3 <- sapply(1:B,function(b){
  if ((b/100==round(b/100)))print(b)
  data3 <- Sterioss[tapply(1:NROW(Sterioss),Sterioss$patnr,function(w){sample(w,1)}),]
  predict(prodlim(Hist(survtime,survstatus)~1,data=data3),times=seq(0,72,12),verbose=F)$surv
})

# Greenwood and Williams
# --------------------------------------------------------------------
X <- prodlim(Hist(survtime,survstatus)~1,data=Sterioss)
x <- prodlim(Hist(survtime,survstatus)~cluster(patnr),data=Sterioss)
y <- predict(x,times=seq(0,72,12))

out <- rbind(c(0,X$se.surv)[y$indices$time+1],
             c(0,x$se.surv)[y$indices$time+1],
             apply(bootx,1,function(x)sqrt(var(x))),
apply(bootx2,1,function(x)sqrt(var(x))),apply(bootx3,1,function(x)sqrt(var(x))))
#apply(bootx2,1,function(x)sd(x)))
dimnames(out) <- list(c("greenwood","cluster","boot.unit","boot.cluster","boot.1unit.cluster"),concat("t",seq(0,72,12)))
out 

apply(bootx,1,function(x)mean(x))

u <- prodlim(Hist(surrogate.time,status)~material+cluster(child),data=primary)

U <- prodlim(Hist(survtime,survstatus)~sex+cluster(patnr),data=Sterioss)

o <- prodlim(Hist(surrogate.time,status)~material+jaw,data=primary)

u <- prodlim(Hist(surrogate.time,status)~material+jaw+cluster(child),data=primary)

v <- clustersurv(Surv(surrogate.time,status)~1|child,
                 data=primary[primary$material=="AM" & primary$jaw=="mandible",])

v <- clustersurv(Surv(surrogate.time,status)~1|child,data=primary)

#comprisk


MVtest$event[MVtest$event==3] <- 0
MVtest <- MVtest[order(MVtest$time,MVtest$event),]

fit1 <- prodlim(Hist(time,event)~1,data=MVtest,subset=group=="B")
cbind(fit1$time,fit1$n.event,fit1$se.cuminc[,1])[fit1$n.event[,1]>0,]

fit2 <- with(MVtest[MVtest$group=="B",],cuminc(time,event))
fit2 <- with(MVtest,cuminc(time,event))
tmp1 <- (fit1$se.cuminc[fit1$n.event[,1]>0,1])^2


tmp2 <- fit2$"1 1"$var[seq(1,length(fit2$"1 1"$var),2)]
cbind(unique(c(0,MVtest$time[MVtest$event==1])),c(0,tmp1),tmp2)




testdat1 <- data.frame(time=c(1,2,2,3,7),status=c(1,2,1,0,2))
z <- prodlim(Hist(time,status)~1,data=testdat1)
Z <- with(testdat1,cuminc(time,status))
library(timereg)
data(TRACE)


a <- prodlim(Hist(time,as.numeric(TRACE$status>0))~1,data=TRACE)
z <- prodlim(Hist(time,status)~1,data=TRACE)



# interval censoring
# --------------------------------------------------------------------

library(Icens)
data(cosmesis)
csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
csub2 <- subset(cosmesis, subset=Trt==1, select=c(L,R))
e1 <- VEM(csub1)
e2 <- VEM(csub2)
par(mfrow=c(2,1))
plot(e1,surv=TRUE)
tmp0 <- PLicens(cosmesis[cosmesis$Trt==0,]$L,cosmesis[cosmesis$Trt==0,]$R,nintervals=0)
lines(tmp0$sieve.L,tmp0$surv,type="s",col=2)
plot(e2,surv=TRUE)
tmp1 <- PLicens(cosmesis[cosmesis$Trt==1,]$L,cosmesis[cosmesis$Trt==1,]$R,nintervals=0)
lines(tmp1$sieve.L,tmp1$surv,type="s",col=2)


# pseudo values
# --------------------------------------------------------------------

set.seed(111)
N <- 1000
testdat <- SimSurv(N=N,surv.dist="rweibull",surv.args=list(shape = 1),surv.baseline=1/100,surv.link="exp",cens.dist="rexp",cens.args=NULL,cens.baseline=1/1000,cens.link="exp",censored=TRUE,cens.max=NULL,keep.uncensored=TRUE,method="transform",verbose=TRUE)
## jackframe <- do.call("rbind",lapply(1:N,function(k){
##   all <- prodlim(Surv(time,status)~1,data=testdat,reverse=FALSE)
##   mink <- prodlim(Surv(time,status)~1,data=testdat[-k,],reverse=FALSE)
##   S <- all$surv
##   Sk <- predictSurv(mink,times=all$time)
##   Jk <- N * S - (N-1) * Sk
##   Jk
## }))


set.seed(111)
N <- 300
testdat <- SimSurv(N=N,surv.dist="rweibull",surv.args=list(shape = 1),surv.baseline=1/100,surv.link="exp",cens.dist="rexp",cens.args=NULL,cens.baseline=1/1000,cens.link="exp",censored=TRUE,cens.max=NULL,keep.uncensored=TRUE,method="transform",verbose=TRUE)
fit <- prodlim(Surv(time,status)~1,data=testdat,reverse=FALSE)
jackframe <- system.time(jackknife.prodlim(fit))
summary((colMeans(jackframe)-fit$surv))


system.time(testpec <- pec(object=fit,
                           formula=Surv(time,status)~1,
                           data=testdat,
                           times=fit$time,
                           replan="boot.632plus",
                           B=300,
                           verbose=FALSE))



N/(N-1)(all$surv-predict(mink,times=all$time)$surv)


set.seed(111)
N <- 10000
k <- 5555
testdat <- SimSurv(N=N,surv.dist="rweibull",surv.args=list(shape = 1),surv.baseline=1/100,surv.link="exp",cens.dist="rexp",cens.args=NULL,cens.baseline=1/1000,cens.link="exp",censored=TRUE,cens.max=NULL,keep.uncensored=TRUE,method=c("transform", "Bender"),verbose=TRUE)
G <- predict(prodlim(Surv(time,status)~1,data=testdat,reverse=TRUE),times=testdat$time,individual=TRUE)
jvk <- N*all$surv-(N-1)*predict(mink,times=all$time)$surv
rhs <- 1-testdat$status[k]*sindex(testdat[k,]$time,testdat$time)/G[k]
plot(unique(testdat$time),jvk,ylim=c(-1,1),type="l")
lines(unique(testdat$time),rhs,type="l",col=2)
abline(v=testdat$time[k])


#set.seed(111)
N <- 10000
testdat <- SimSurv(N=N,surv.dist="rweibull",surv.args=list(shape = 1),surv.baseline=1/100,surv.link="exp",cens.dist="rexp",cens.args=NULL,cens.baseline=1/1000,cens.link="exp",censored=F,cens.max=NULL,keep.uncensored=TRUE,method=c("transform", "Bender"),verbose=TRUE)
k <- 5555
all <- prodlim(Surv(time,status)~1,data=testdat,reverse=F)
tk <- testdat[-k,]
mink <- prodlim(Surv(time,status)~1,data=tk,reverse=F)
print(max(all$surv-predict(mink,times=all$time)$surv))
pvk <- (N-1)*(all$surv-predict(mink,times=all$time)$surv)
plot(all$time,pvk,ylim=c(-1,1),xlim=c(0,10))
abline(v=testdat$time[k])



# plot.Hist
# --------------------------------------------------------------------

## two-state model

frame2 <- data.frame(time=c(1,2),status=c(1,1))
with(frame2,plot(Hist(time,status)))
