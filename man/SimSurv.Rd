\name{SimSurv}
\alias{SimSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating survival data }
\description{
  Censored event times are drawn from user specified conditional distributions
  given simulated covariates.
}
\usage{
SimSurv(N,
	surv,
	cens,
	cova,
	verbose=1,
	...)
}
\arguments{
  \item{N}{Sample size}
  \item{surv}{Dummy argument. The survival distribution is determined by
    arguments of the form \code{surv.arg}. See Details}
  \item{cens}{If FALSE data are left uncensored. Otherwise the censoring
    distribution is specified by using arguments of the form
    \code{cens.arg}. See Details.}
  \item{cova}{either a matrix with \code{N} rows,
    or a named list to generate the covariates. each entry
    is a list where the first element is the names of the
    function used to draw the covariate values and the remaining
    elements are arguments passed to that function. For example
    \code{cova=list(X = list(dist = "rlnorm", meanlog = 2, sdlog = 0.4))}
    generates a single log-normally distributed covariate called
    \code{X}.
  }
  \item{verbose}{Set to \code{FALSE} to shut up in simulations}
  \item{...}{used for convenient argument specification, e.g.
    \code{surv.transform.X2=function(x)x^2}
    will overwrite a corresponding entry
    of \code{surv} for the \code{transform} of the
    covariate \code{X2}
  }
}
\details{
  \code{surv.dist} name of the function to draw the survival distribution
  \code{surv.args} list of extra arguments to \code{surv.dist}.
  \code{surv.baseline} the baseline risk 
  \code{surv.link} names of the link to the covariates
  \code{surv.coef} numeric vector of regression coefficients. the ith is
    used for the ith entry of \code{cova}
  \code{surv.transform} list of function names one for each covariate. the ith is applied
    to values of the ith covariate before it is linked to generate the
    survival times.
  \code{cens.dist} name of the function to draw the censoring distribution
  \code{cens.args} list of extra arguments to \code{dist}
  \code{cens.baseline} the baseline risk
  \code{cens.link} names of the link to the covariates
  \code{cens.max} maximal value where all event times are right censored
  \code{cens.type} "right" for right censored data, "interval" for interval
    censored data
  \code{cens.coef} numeric vector of regression coefficients. the ith is
    used for the ith entry of \code{cova}
  \code{cens.transform}list of function names one for each covariate,
    where the ith element is applied
  to values of the ith covariate before it is linked to generate the
  survival times.
  
  Possible distributions to generate covariates:

  X.lognorm = list(dist = "rlnorm", meanlog = 2, sdlog = 0.4)
   
   X.unif = list(dist = "runif", min = 0, max = 10),
   
   X.exp = list(dist = "rexp", rate = 0.4)
   
   X.bernoulli = list(dist = "rbinom", size = 1, prob = 0.3)
   
   X.binom = list(dist = "rbinom", size = 4, prob = 0.8)
   
   X.nbinom = list(dist = "rnbinom", size = 3, mu = 2)
   
   X.poisson = list(dist = "rpois", lambda = 1.3) 
}
\value{
  A data.frame:
  \item{time}{the right censored event times}
  \item{status}{the survival status}
  \item{X}{the values of the covariate \code{X}}
  \item{f.X}{the transformed values of the covariate \code{X}}
  \item{time}{the uncensored event times}

  Attributes of the data.frame
  \code{formula} a formula to evaluate the generated event history object
}
\references{

  Ralf Bender, Thomas Augustin, and Maria Blettner. Generating 
survival times to simulate Cox proportional hazards models by Ralf Bender, 
Thomas Augustin and Maria Blettner, Statistics in Medicine 2005; 
24:1713-1723. Stat Med, 25(11):1978-9, 2006.

}
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
% \note{
% }
\seealso{\code{\link{prodlim}}}
\examples{
\donttest{
SimSurv(10)
SurvData=SimSurv(100,cens.baseline=1/10,surv.baseline=2)
Hist(SurvData$time,SurvData$status)
prodlim(Hist(time,status)~1,data=SurvData)
plot(prodlim(Hist(time,status)~1,data=SurvData))

SurvData=SimSurv(100,cens.baseline=1/10,surv.baseline=2,surv.coef=c(-1,-2),
         cova=list( X.exp = list(dist = "rexp", rate = 0.4),
         X.bernoulli = list(dist = "rbinom", size = 1, prob = 0.3)))


SurvData=SimSurv(100,cens.baseline=1/10,surv.baseline=2,surv.coef=c(-1,-2),cens.coef=c(0,1),
          cova=list( X.exp = list(dist = "rexp", rate = 0.4),
          X.bernoulli = list(dist = "rbinom", size = 1, prob = 0.3)))

set.seed(8)
SurvFrame <- SimSurv(N=100,
		     surv=list(dist = "rweibull",args = list(shape=1),baseline = 1/100,link="exp",coef=1),
		     cens=list(dist = "rexp",args=NULL,baseline=1/1000,link="exp",coef=0,max = 150,type = "interval",lateness=1,unit=10),
		     cova=list(Age = list("rnorm", mean = 0, sd = 2),Sex = list("rbinom", size=1,prob=.5)),
		     keep.uncensored=TRUE,
		     method=c("simulation"),
		     verbose=1)

library(survival)
coxph(Surv(time,status==0)~X.exp+X.bernoulli,data=SurvData)

coxph(Surv(time,status)~X.exp+X.bernoulli,data=SurvData)

# Simulate survival times based on a parametric survival fit

library(survival)
data(ovarian)
fitSurv <- survreg(Surv(futime,fustat)~age+rx,data=ovarian)
fitCens <- survreg(Surv(futime,1-fustat)~age+rx,data=ovarian)
# note: survreg fits the accelerated failure time model version of the
#       Weibull survival regression model. 
#       The parameters need to be transformed in order to
#       simulate from a fitted model using SimSurv:
survCoef <- -fitSurv$coef[-1]/fitSurv$scale
survBaseline <- exp(-fitSurv$coef["(Intercept)"]/fitSurv$scale)
survShape <- 1/fitSurv$scale

censCoef <- -fitCens$coef[-1]/fitCens$scale
censBaseline <- exp(-fitCens$coef["(Intercept)"]/fitCens$scale)
censShape <- 1/fitCens$scale
X <- model.frame(~age+rx,data=ovarian)

simuData <- SimSurv(N=26,
		    cova=X,
		    surv.coef=survCoef,
		    surv.baseline=survBaseline,
		    surv.shape=survShape,
                    cens.model="Cox-Weibull",
		    cens.coef=censCoef,
		    cens.baseline=censBaseline,
		    cens.shape=censShape)

simuFit <- survreg(Surv(time,status)~age+rx,data=simuData)
cbind(coef(simuFit),coef(fit))

plot(prodlim(Surv(time,status)~rx,data=simuData),confint=F)
plot(prodlim(Surv(futime,fustat)~rx,data=ovarian),confint=F,add=TRUE,lty=2)

## another example

surv.coef=c(0.1, 0.03, -0.5, 0.5) 
cens.coef=c(0.1, 0.03,-0.5, 0.5) 
variables=list( 
     x1=list(dist="rnorm", mean=26, sd=2), 
     x2=list(dist="rnorm", mean=50, sd=7), 
     x3=list(dist="rbinom", size=1, prob=0.2), 
     x4=list(dist="rbinom", size=1, prob=0.5)) 
surv.base<-0.00001 
cens.base<-0.00009 
sim.dat<-SimSurv(N=500,
		 cova=variables,
		 surv.coef=surv.coef,
		 cens.coef=cens.coef,
		 surv.model="Cox-exponential",
		 cens.model="Cox-exponential",
		 surv.baseline=surv.base,
		 cens.baseline=cens.base,
		 verbose=F)

## adjusting the censoring percentage
## for a given survival setting

SNP1 <- rbinom(181,2,.3)
stage <- rbinom(181,2,.5)
setup  <- list(N=181,surv.baseline=1/100, cova=list(a=SNP1,b=stage), surv.coef=c(1.3,1.8), cens.baseline=1/100)
d <- do.call("SimSurv",setup)
set30 <- prodlim:::find.baseline(x=.3,setup)
d30 <- do.call("SimSurv",set30)
set70 <- prodlim:::find.baseline(x=.7,setup)
d70 <- do.call("SimSurv",set70)

## to visulize the shape of the censoring 
## distribution do this:

plot(prodlim(Hist(time,status)~1,data=d30,rev=TRUE))
plot(prodlim(Hist(time,status)~1,data=d70,rev=TRUE))

## what matters is (not the mean censoring time but)
## the support of the censoring distribution,
## that is the time region where the probability of
## being uncensored is different from zero


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% at least one, from doc/KEYWORDS
