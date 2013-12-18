library(lava)
library(prodlim)
check.code(prodlim)

## setting 1: predictor Xp is associated to
##            randomly chosen np other predictors
##            the effect sizes are a random numbers
P <- 100
sparseness <- 8
covnames <- paste("X",1:P,sep="")
## 1: generat the covariate matrix
set.seed(13234)
covform <- lapply(1:(P-1),function(p){
  ## find predictors correlated with Xp 
  selectX <- paste("X",sample(1:p,size=min(rbinom(1,sparseness,.4),p),replace=FALSE),sep="")
  if ((length(selectX)==1) && selectX=="X"){## no correlation
    formula(paste(covnames[p+1],"~f(X1,0)"))
  }
  else{
    ## find coefficients for the association
    coefX <- round(rnorm(length(selectX),mean=0,sd=0.5),2)
    formula(paste(covnames[p+1],"~",paste(paste("f(",selectX,",",coefX,")"),collapse="+")))
  }
})
## the setting up is a bit slow when P gets larger:
model <- lvm(covform)

## the resulting graphical model:
plot(model)
## the generated correlation structure:
model$fix
## this works as expected:
if (FALSE){
  d <- sim(model,1000)
  testform <- formula(paste(paste("X",P,"~",sep=""),paste(paste("X",c(1:(P-1)),sep=""),collapse="+")))
  summary(glm(testform,data=d))
}
## now competing risks data
## variables X1, X12, X13, X17, X22, X48 have a direct effect on
## competing cause 1
## variables X17, X22, X48, X81, X82, X99 have a direct effect on
## competing cause 2
## hence: variables X17, X22, X48 have a direct effect on
## both competing causes
compriskEffects <- list(X1=list(cr1=-.3,cr2=0),
                        X12=list(cr1=.7,cr2=0),
                        X13=list(cr1=-.5,cr2=0),
                        X17=list(cr1=.7,cr2=.4),
                        X22=list(cr1=-.6,cr2=-.4),
                        X48=list(cr1=-.5,cr2=.4),
                        X81=list(cr1=0,cr2=.4),
                        X82=list(cr1=0,cr2=-.7),
                        X99=list(cr1=0,cr2=.9))
cause1Effects <- sapply(compriskEffects,function(x){x[[1]]})
cause1Effects <- cause1Effects[cause1Effects!=0]
cause2Effects <- sapply(compriskEffects,function(x){x[[2]]})
cause2Effects <- cause2Effects[cause2Effects!=0]
cause1Formula <- formula(paste("Time1~",paste(paste("f(",names(cause1Effects),",",cause1Effects,")",sep=""),collapse="+")))
cause2Formula <- formula(paste("Time2~",paste(paste("f(",names(cause2Effects),",",cause2Effects,")",sep=""),collapse="+")))
regression(model) <- cause1Formula
regression(model) <- cause2Formula
distribution(model,~Time1) <- weibull.lvm()
distribution(model,~Time2) <- weibull.lvm()
## show direct and indirect effects:
plot(model)

## simulate competing risks data
set1CompRiskData <- function(N){
  data <- sim(model,N)
  data$time <- pmin(data$Time1,data$Time2)
  data$event <- 1
  data$event[data$Time2<data$Time1] <- 2
  data
}
testForm1 <- formula(paste("Surv(time,event==1)~",paste(paste("X",1:P,sep=""),collapse="+")))
testForm2 <- formula(paste("Surv(time,event==2)~",paste(paste("X",1:P,sep=""),collapse="+")))

## low dimensional setting (N=1000, P=100)
d <- set1CompRiskData(1000)
## full Cox models
coxph(testForm1,data=d)
coxph(testForm2,data=d)
## backward elimination
library(rms)
fastbw(cph(testForm1,data=d))
fastbw(cph(testForm2,data=d))

## low dimensional setting (N=1000, P=100)
set.seed(232)
d <- set1CompRiskData(100)
## Cox models cannot fit
coxph(testForm1,data=d)
coxph(testForm2,data=d)
library(randomForestSRC)
cause1Forest <- rfsrc(testForm1,data=d)



