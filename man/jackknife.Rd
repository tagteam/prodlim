\name{jackknife}
\alias{jackknife}
\alias{jackknife.competing.risks}
\alias{jackknife.survival}
\title{Compute jackknife pseudo values.}
\usage{
jackknife(object, times, cause, keepResponse = FALSE, ...)
}
\arguments{
  \item{object}{Object of class \code{"prodlim"}.}

  \item{times}{Time points at which to compute pseudo
  values.}

  \item{cause}{For competing risks the cause of failure.}

  \item{keepResponse}{If \code{TRUE} add the model
  response, i.e. event time, event status, etc. to the
  result.}

  \item{...}{not used}
}
\description{
Compute jackknife pseudo values.
}
\details{
Compute jackknife pseudo values based on marginal
Kaplan-Meier estimate of survival, or based on marginal
Aalen-Johannsen estimate of cumulative incidence.
}
\note{
The R-package pseudo does a similar job, and appears to be
a little faster in small samples, but much slower in large
samples. See examples.
}
\examples{
## pseudo-values for survival models

d=SimSurv(20)
f=prodlim(Hist(time,status)~1,data=d)
jackknife(f,times=c(30,50))

## in some situations it may be useful to attach the
## the event time history to the result
jackknife(f,times=c(30,50),keepResponse=TRUE)

# pseudo-values for competing risk models
d=SimCompRisk(10)
f=prodlim(Hist(time,event)~1,data=d)
jackknife(f,times=c(3,10),cause=1)
jackknife(f,times=c(3,10,17),cause=2)

\donttest{
library(prodlim)
library(pseudo)

# comparison to pseudoci
# make sure we get the same
# results with both packages

set.seed(17)
N <- 200
ddd <- SimCompRisk(200)
ttt <- c(3,5,10)
# ttt <- ddd$time
fff <- prodlim(Hist(time,event)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})

## check individual 2
all(round(jack2[,2],9)==round(jack[[3]]$cause1[,2],9))

## check all individuals
all(sapply(1:N,function(x){
a <- round(jack[[3]]$cause1[x,],8)
b <- round(jack2[x,],8)
# all(a[!is.na(a)]==b[!is.na(b)])
all(a[!is.na(a)]==b[!is.na(a)])
}))


## the pseudoci function seems only slightly slower
## for small sample sizes (up to ca. 200) but
## much slower for large sample sizes:

set.seed(17)
N <- 200
ddd <- SimCompRisk(200)
ttt <- c(3,5,10)
# ttt <- ddd$time
fff <- prodlim(Hist(time,event)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,event,ttt))[,-c(1:2,seq(4,N+2,2))])
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
all(round(jack2[,1],9)==round(jack[,1],9))

set.seed(17)
N <- 2000
ddd <- SimCompRisk(2000)
ttt <- c(3,5,10)
fff <- prodlim(Hist(time,event)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
}
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\references{
Andersen PK & Perme MP (2010). Pseudo-observations in
survival analysis Statistical Methods in Medical Research,
19(1), 71-99.
}
\seealso{
\code{\link{prodlim}}
}
\keyword{survival}

