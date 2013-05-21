## The function `Hist' (short for history) generalizes
## the function `Surv' which is defined in the survival package

# {{{ two state survival model 

## If the status variable takes only two values, then `Surv' and
## `Hist' do more or less the same. The attributes and the
## print/summary functions are different. 

library(prodlim)
library(survival)
set.seed(17)
dat <- SimSurv(100)
sH1 <- Hist(dat$time,dat$status)
sH2 <- Surv(dat$time,dat$status)
stopifnot(all(unclass(sH1)==unclass(sH2)))

## summary(sH1)
## summary(sH2)
## without(formula)
a <- 1:10
b <- rep(0,5)
prodlim(Surv(a,b)~1)

set.seed(10)
dat <- SimSurv(100)
d <- SimSurv(100)
f <- prodlim(Hist(time,status)~1,data=d)
sum <- summary(f,times=d$time)

library(survival)
fcontrol <- survfit(Surv(time,status)~1,data=d)
sumcontrol <- summary(fcontrol,times=d$time)

# second test: a competing risk model 
# --------------------------------------------------------------------

## require(prodlim)
## require(cmprsk)

## x <- 1:15
## y <- c(2,1,0,1,2,1,0,1,1,0,1,1,2,1,1)
## A <- cuminc(x,y)

## Hist(x,y)
## B <- prodlim(Hist(x,y)~1)
#plot(B)
#plot(A)
#plot(B,cause=c(1,2))

## all(A$"1 1"$est == predict(B,times=c(0,A$"1 1"$time[-length(A$"1 1"$time)]))[[1]])
## all(A$"1 2"$est == predict(B,cause=2,times=c(0,A$"1 2"$time[-length(A$"1 2"$time)]))[[1]])

## test <- data.frame(time=rep(1:10,2),status=rep(0:1,10))
## pfit <- prodlim(Hist(time,status)~1,data=test)
## summary(pfit,times=c(2,4.4,9))


## prodlim(Surv(time,status)~1,data=test) ## same as pfit
## fit <- survfit(Surv(time,status)~1,data=test) ## same as pfit$surv
## sfit <- summary(fit,times=c(2,4.4,9))
## all(round(pfit$surv,10)==round(survfit(Surv(time,status)~1,data=test)$surv,10))

# survival
## require(survival)
pfit <- prodlim(Hist(time,status)~1,data=dat)
all(round(pfit$surv,10)==round(survfit(Surv(time,status)~1,data=dat)$surv,10))

dat <- SimSurv(100)
a <- prodlim(Hist(time,status)~1,data=dat)
## b <- prodlim(Surv(time,status)~1,data=dat,reverse=TRUE)
## c <- prodlim(Surv(time,status)~X2,data=dat,reverse=TRUE)
d <- prodlim(Surv(time,status)~NN(X1),data=dat,reverse=TRUE)

## if (!is.function("cluster")) cluster <- function(x)x
## class(c) <- "dynamic"
#cluster
## testdat <- data.frame(time=c(1,2,2,3,7),status=c(0,1,1,0,1),patnr=c(1,2,1,3,2))
## y <- prodlim(Hist(time,status)~cluster(patnr),data=testdat)

## library(dental)
## x <- clustersurv(Surv(time,status)~1|patnr,data=testdat)
## x <- prodlim(Hist(survtime,survstatus)~cluster(patnr),data=Sterioss)
## data(primary)
## data(Sterioss)

## w <- as.numeric(Sterioss$patnr)


# bootstrapping study units
# --------------------------------------------------------------------
## N <- NROW(Sterioss)
## boot <- matrix(sapply(1:B,function(b){sample(1:N,replace=TRUE)}),nrow=N,ncol=B)
## bootx <- sapply(1:B,function(b){
  ## if ((b/100==round(b/100)))print(b)
  ## u <- predict(prodlim(Hist(survtime,survstatus)~1,data=Sterioss[boot[,b],]),times=seq(0,72,12),verbose=F)$surv
## })
# bootstrapping patients
# --------------------------------------------------------------------
## B <- 1000
## N2 <- length(unique(Sterioss$patnr))
## Sterioss.patlist <- split(Sterioss[,c("survstatus","survtime","patnr")],Sterioss$patnr)
## boot2 <- matrix(sapply(1:B,function(b){sample(1:N2,replace=TRUE)}),nrow=N2,ncol=B)
## bootx2 <- sapply(1:B,function(b){
  ## if ((b/100==round(b/100)))print(b)
  ## data2 <- do.call("rbind",Sterioss.patlist[boot2[,b]])
  ## if(sum(data2$survstatus)==0) print(paste(b,"failed"))
  ## else predict(prodlim(Hist(survtime,survstatus)~1,data=data2),times=seq(0,72,12),verbose=F)$surv
## })
## fit <- prodlim(Surv(survtime,survstatus)~cluster(patnr),data=Sterioss)
## test <- bootCluster.product.limit(fit,
                                  ## times=sort(unique(Sterioss$survtime)),
                                  ## B=10)


# repeated sampling one unit per patient 
# --------------------------------------------------------------------
## B <- 1000
## bootx3 <- sapply(1:B,function(b){
  ## if ((b/100==round(b/100)))print(b)
  ## data3 <- Sterioss[tapply(1:NROW(Sterioss),Sterioss$patnr,function(w){sample(w,1)}),]
  ## predict(prodlim(Hist(survtime,survstatus)~1,data=data3),times=seq(0,72,12),verbose=F)$surv
## })

# Greenwood and Williams
# --------------------------------------------------------------------
## X <- prodlim(Hist(survtime,survstatus)~1,data=Sterioss)
## x <- prodlim(Hist(survtime,survstatus)~cluster(patnr),data=Sterioss)
## y <- predict(x,times=seq(0,72,12))

## out <- rbind(c(0,X$se.surv)[y$indices$time+1],
             ## c(0,x$se.surv)[y$indices$time+1],
             ## apply(bootx,1,function(x)sqrt(var(x))),
## apply(bootx2,1,function(x)sqrt(var(x))),apply(bootx3,1,function(x)sqrt(var(x))))
## #apply(bootx2,1,function(x)sd(x)))
## dimnames(out) <- list(c("greenwood","cluster","boot.unit","boot.cluster","boot.1unit.cluster"),concat("t",seq(0,72,12)))
## out 

## apply(bootx,1,function(x)mean(x))

## u <- prodlim(Hist(surrogate.time,status)~material+cluster(child),data=primary)
## U <- prodlim(Hist(survtime,survstatus)~sex+cluster(patnr),data=Sterioss)
## o <- prodlim(Hist(surrogate.time,status)~material+jaw,data=primary)
## u <- prodlim(Hist(surrogate.time,status)~material+jaw+cluster(child),data=primary)
## v <- clustersurv(Surv(surrogate.time,status)~1|child,
## data=primary[primary$material=="AM" & primary$jaw=="mandible",])
## v <- clustersurv(Surv(surrogate.time,status)~1|child,data=primary)


# interval censoring
# --------------------------------------------------------------------

## library(Icens)
## data(cosmesis)
## csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
## csub2 <- subset(cosmesis, subset=Trt==1, select=c(L,R))
## e1 <- VEM(csub1)
## e2 <- VEM(csub2)
## par(mfrow=c(2,1))
## plot(e1,surv=TRUE)
## tmp0 <- PLicens(cosmesis[cosmesis$Trt==0,]$L,cosmesis[cosmesis$Trt==0,]$R,nintervals=0)
## lines(tmp0$sieve.L,tmp0$surv,type="s",col=2)
## plot(e2,surv=TRUE)
## tmp1 <- PLicens(cosmesis[cosmesis$Trt==1,]$L,cosmesis[cosmesis$Trt==1,]$R,nintervals=0)
## lines(tmp1$sieve.L,tmp1$surv,type="s",col=2)










