\name{jackknife}
\alias{jackknife}
\alias{leaveOneOut}
\alias{leaveOneOut.survival}
\alias{leaveOneOut.competing.risks}
\alias{jackknife.survival}
\alias{jackknife.competing.risks}
\title{
Compute jackknife pseudo values.
}
\description{
Compute jackknife pseudo values based on marginal Kaplan-Meier estimate
of survival, or based on marginal Aalen-Johannsen estimate of cumulative incidence.
}
\usage{
jackknife(object, times,keepResponse=FALSE, ...)
jackknife.survival(object, times,keepResponse=FALSE, ...)
jackknife.competing.risks(object, times, cause, keepResponse=FALSE, ...)
leaveOneOut(object, times, ...)
leaveOneOut.survival(object, times, lag=FALSE, ...)
leaveOneOut.competing.risks(object, times, cause, ...)
}
\arguments{
  \item{object}{
    Object of class \code{"prodlim"}.
  }
  \item{times}{
    Time points at which to compute pseudo values.
  }
  \item{cause}{
    For competing risks the cause of failure.
  }
  \item{keepResponse}{
    If \code{TRUE} add the model response, i.e. event time, event
    status, etc. to the result.
  }
  \item{lag}{
    Logical. If \code{TRUE} lag the result, i.e. compute S(t-) instead
    of S(t). 
  }
  \item{\dots}{
  }
}
% \details{
% }
% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
\references{
Andersen PK & Perme MP (2010). Pseudo-observations in survival analysis
Statistical Methods in Medical Research, 19(1), 71-99.


}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
 \note{
The R-package pseudo does a similar job, and appears to be a little faster in
small samples, but much slower in large samples. See examples.
 }

\seealso{
  \code{\link{prodlim}}
}
\examples{
\donttest{

## pseudo-values for survival models

d=SimSurv(100) 
f=prodlim(Hist(time,status)~1,data=d) 
jackknife(f,times=c(30,50))

## in some situations it may be useful to attach the
## the event time history to the result
jackknife(f,times=c(30,50),keepResponse=TRUE)

# pseudo-values for competing risk models
d=prodlim:::SimCompRisk(100) 
f=prodlim(Hist(time,cause)~1,data=d) 
jackknife(f,times=c(3,10),cause=1)
jackknife(f,times=c(3,10,17),cause=2)

# comparison to pseudoci
# make sure we get the same
# results with both packages

library(prodlim)
library(pseudo)

set.seed(17)
N <- 200
ddd <- data.frame(time=1:N,cause=rbinom(N,2,.5),X=rbinom(N,1,.5))
ttt <- c(3,5,10)
# ttt <- ddd$time
fff <- prodlim(Hist(time,cause)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,cause,ttt))[,-c(1:2,seq(4,N+2,2))])
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})

## check individual 1
all(round(jack2[,1],9)==round(jack[,1],9))

## check all individuals
all(sapply(1:N,function(x){
a <- round(jack[x,],8)
b <- round(jack2[x,],8)
# all(a[!is.na(a)]==b[!is.na(b)])
all(a[!is.na(a)]==b[!is.na(a)])
}))


## the pseudoci function seems only slightly slower
## for small sample sizes (up to ca. 200) but
## much slower for large sample sizes:

set.seed(17)
N <- 200
ddd <- data.frame(time=1:N,cause=rbinom(N,2,.5),X=rbinom(N,1,.5))
ttt <- c(3,5,10)
# ttt <- ddd$time
fff <- prodlim(Hist(time,cause)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,cause,ttt))[,-c(1:2,seq(4,N+2,2))])
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
all(round(jack2[,1],9)==round(jack[,1],9))

set.seed(17)
N <- 2000
ddd <- data.frame(time=1:N,cause=rbinom(N,2,.5),X=rbinom(N,1,.5))
ttt <- c(3,5,10)
# ttt <- ddd$time
fff <- prodlim(Hist(time,cause)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,cause,ttt))[,-c(1:2,seq(4,N+2,2))])
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
all(round(jack2[,1],9)==round(jack[,1],9))

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
