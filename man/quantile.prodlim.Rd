\name{quantile.prodlim}
\alias{quantile.prodlim}
\title{Quantiles for Kaplan-Meier and Aalen-Johansen estimates.}
\usage{
\method{quantile}{prodlim}(x, q, cause = 1, ...)
}
\arguments{
  \item{x}{Object of class \code{"prodlim"}.}

  \item{q}{Quantiles. Vector of values between 0 and 1.}

  \item{cause}{For competing risks the cause of interest.}

  \item{...}{not used}
}
\description{
Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
}
\examples{
set.seed(1)
d=SimSurv(30)
f=prodlim(Hist(time,status)~X1,data=d)
# median
quantile(f,.5)
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\keyword{survival}

