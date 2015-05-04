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
Aalen-Johansen estimate of cumulative incidence.
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
jackknife(f,times=c(3,5))

## in some situations it may be useful to attach the
## the event time history to the result
jackknife(f,times=c(3,5),keepResponse=TRUE)

# pseudo-values for competing risk models
d=SimCompRisk(10)
f=prodlim(Hist(time,event)~1,data=d)
jackknife(f,times=c(3,10),cause=1)
jackknife(f,times=c(3,10,17),cause=2)
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

