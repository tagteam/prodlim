\name{prodlim}
\alias{prodlim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{product limit method }
\description{
  Nonparametric estimation in event history analysis. Featuring fast
  algorithms and user friendly syntax adapted from the survival package.
  The product limit algorithm is used for right censored data; the
  self-consistency algorithm for interval censored data.
}

\usage{prodlim(formula,
	       data=parent.frame(),
	       subset,
	       na.action,
	       reverse=FALSE,
	       conf.int=0.95,
	       bandwidth=NULL,
	       discrete.level=3,
               maxiter=1000,
               grid,
               tol=7,
               method=c("npmle","one.step","impute.midpoint","impute.right"),
               exact=TRUE)
     }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula whose left hand side is 
                 a \code{Hist} object. In some special cases it can also be a
                 \code{Surv} response object, see the details section. The right hand side is as
                 usual a linear combination of covariates which may
                 contain at most one continuous factor. Whether or not a covariate is
                 recognized as continuous or discrete depends on its class and on the
                 argument \code{discrete.level}. The right hand side may also be used
                 to specify clusters, see the details section.
               }
  \item{data}{
    A data.frame in which all the variables of
    \code{formula} can be interpreted.
  }
  \item{subset}{
    Expression identifying a subset of the rows of the
    data for the fitting process.
  }
  \item{na.action}{
    A missing-data filter function, applied to the model.frame,
    after any subset argument has been used.  Default is
    \code{options()$na.action}.
  }
  \item{reverse}{
    For right censored data, if reverse=TRUE then the censoring
    distribution is estimated.
  }
  \item{conf.int}{
    The level (between 0 and 1) for two-sided pointwise confidence
    intervals. Defaults to 0.95.
  }
  \item{bandwidth}{
    Smoothing parameter for symmetric nearest neighborhoods
    based on the values of a continuous covariate. See function
    \code{neighborhood} for details.
  }
  \item{discrete.level}{
    Numeric covariates are treated as factors when their number
    of unique values exceeds not \code{discrete.level}. Otherwise
    the product limit method is applied, in overlapping neighborhoods
    according to the bandwidth.
  }
  \item{maxiter}{
    For interval censored data only.
    Maximal number of iterations to obtain the nonparametric maximum likelihood estimate.
    Defaults to 1000.
  }
  \item{grid}{
    For interval censored data only.
    When method=one.step grid for one-step product limit estimate. Defaults to
    sorted list of unique left and right endpoints of the observed intervals.
  }
  \item{tol}{
    For interval censored data only.
    Numeric value whose negative exponential is used as convergence criterion
    for finding the nonparametric maximum likelihood estimate.
    Defaults to 7 meaning exp(-7).}
  \item{method}{
    For interval censored data only.
    If equal to \code{"npmle"} (the default) use the usual Turnbull algorithm, else
    the product limit version of the
    self-consistent estimate.
  }
  \item{exact}{If TRUE the grid of time points used for estimation
    includes all the L and R endpoints of the observed intervals.}
}
\details{

  The response of \code{formula} (ie the left hand side of the `~'
  operator) specifies the model.
  
  In two-state models -- the classical survival case -- the standard
  Kaplan-Meier method is applied.  For this the response can be
  specified as a \code{\link{Surv}} or as a \code{\link{Hist}}
  object. The \code{\link{Hist}} function allows you to change the code
  for censored observations,
  e.g. \code{Hist(time,status,cens.code="4")}. 
  
  Besides a slight gain of computing efficiency, there are some extensions that
  are not included in the current version of the survival package:

  (0) The Kaplan-Meier estimator for the censoring times
  \code{reverse=TRUE} is correctly
  estimated when there are ties between event and censoring
  times.

  (1) A conditional version of the Kaplan-Meier estimator for at most
  one continuous predictors using symmetrical nearest neighborhoods
  (Beran 1981, Stute 1984, Akritas 1994).

  (2) For cluster-correlated data the right hand side of
  \code{formula} may identify a \code{\link{cluster}} variable. In that case
  Greenwood's variance formula is replaced by the formula of Ying \& Wei
  (1994).

  (3) Competing risk models can be specified via
  \code{\link{Hist}} response objects in \code{formula}.

  The Aalen-Johannsen estimator is applied for estimating the cumulative
  incidence functions for all causes.
  The advantage over the function \code{cuminc} of the cmprsk
  package are user-friendly model specification via \code{\link{Hist}}
  and sophisticated print, summary, predict and plot methods.

  Under construction:

  (U0) Interval censored event times specified via
  \code{\link{Hist}} are used to find the nonparametric maximum
  likelihood estimate. Currently this works only for two-state models
  and the results should match with those from the package `Icens'.
  
  (U1) Extensions to more complex multi-states models
  
  (U2) The nonparametric maximum likelihood estimate for interval
       censored observations of competing risks models.
  
}
\value{
  Object of class "prodlim" for which there are print, predict, life.table, summary and plot
  methods.
}
\references{Andersen, Borgan, Gill, Keiding (1993) Springer
  `Statistical Models Based on Counting Processes'

  Akritas (1994) The Annals of Statistics 22, 1299-1327
 Nearest neighbor estimation of a bivariate distribution under random censoring.

  R Beran (1981) http://anson.ucdavis.edu/~beran/paper.html
  `Nonparametric regression with randomly censored survival data' 

  Stute (1984) The Annals of Statistics 12, 917--926
  `Asymptotic Normality of Nearest Neighbor Regression Function
  Estimates'

  Ying, Wei (1994) Journal of Multivariate Analysis 50, 17-29
  The Kaplan-Meier estimate for dependent failure time observations
  
}

\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}

%  ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{predictSurv}},
  \code{\link{predictSurvIndividual}},
  \code{\link{predictCuminc}},
  \code{\link{Hist}},
  \code{\link{neighborhood}},
  \code{\link{Surv}},
  \code{\link{survfit}},
  \code{\link{strata}},
  \code{\link{cluster}},
  \code{\link{NN}}}
\examples{

##---------------------two-state survival model------------

dat <- SimSurv(100)
with(dat,plot(Hist(time,status)))
fit <- prodlim(Hist(time,status)~1,data=dat)
print(fit)
plot(fit)
summary(fit)

## --------------------clustered data---------------------
library(survival)
dat <- data.frame(time=rexp(100),status=rbinom(100,1,.3),X=rbinom(100,1,.5),Z=rnorm(100,10,3),patnr=sample(1:10,size=100,replace=TRUE))
fit <- prodlim(Hist(time,status)~cluster(patnr),data=dat)
print(fit)
plot(fit)
summary(fit)


##-----------compare Kaplan-Meier to survival package---------

dat <- SimSurv(100)
pfit <- prodlim(Surv(time,status)~1,data=dat)
pfit <- prodlim(Hist(time,status)~1,data=dat) ## same thing
sfit <- survfit(Surv(time,status)~1,data=dat,conf.type="plain")
##  same result for the survival distribution function 
all(round(pfit$surv,12)==round(sfit$surv,12))
summary(pfit,digits=3)
summary(sfit,times=quantile(unique(dat$time)))

##------------------censoring survival function----------------

rpfit <- prodlim(Hist(time,status)~1,data=dat,reverse=TRUE)
rsfit <- survfit(Surv(time,1-status)~1,data=dat,conf.type="plain") 
## not the same if there are ties!
all(round(rpfit$surv,3)==round(rsfit$surv,3))

##-------------------stratified Kaplan-Meier---------------------

pfit.X2 <- prodlim(Surv(time,status)~X2,data=dat)
summary(pfit.X2)
summary(pfit.X2,intervals=TRUE)
plot(pfit.X2)

##----------continuous covariate: Stone-Beran estimate------------

prodlim(Surv(time,status)~X1,data=dat)

##-------------both discrete and continuous covariates------------

prodlim(Surv(time,status)~X2+X1,data=dat)

##----------------------interval censored data----------------------

dat <- data.frame(L=1:10,R=c(2,3,12,8,9,10,7,12,12,12),status=c(1,1,0,1,1,1,1,0,0,0))
with(dat,Hist(time=list(L,R),event=status))

dat$event=1
npmle.fitml <- prodlim(Hist(time=list(L,R),event)~1,data=dat)

##-------------competing risks-------------------

CompRiskFrame <- data.frame(time=1:100,event=rbinom(100,2,.5),X=rbinom(100,1,.5))
crFit <- prodlim(Hist(time,event)~X,data=CompRiskFrame)
summary(crFit)
plot(crFit)
summary(crFit,cause=2)
plot(crFit,cause=2)


# Changing the cens.code:
dat <- data.frame(time=1:10,status=c(1,2,1,2,5,5,1,1,2,2))
fit <- prodlim(Hist(time,status)~1,data=dat)
print(fit$model.response)
fit <- prodlim(Hist(time,status,cens.code="2")~1,data=dat)
print(fit$model.response)
plot(fit)
plot(fit,cause="5")
}

\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{nonparametric}
\keyword{cluster}% 
