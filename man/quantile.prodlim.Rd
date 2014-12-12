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
library(lava)
set.seed(1)
d=SimSurv(30)
f=prodlim(Hist(time,status)~1,data=d)
f1=prodlim(Hist(time,status)~X1,data=d)
# default: median and IQR
quantile(f)
quantile(f1)
# median alone
quantile(f,.5)
quantile(f1,.5)

# competing risks
set.seed(3)
dd = SimCompRisk(30)
ff=prodlim(Hist(time,event)~1,data=dd)
ff1=prodlim(Hist(time,event)~X1,data=dd)
## default: median and IQR
quantile(ff)
quantile(ff1)

print(quantile(ff1),na.val="NA")
print(quantile(ff1),na.val="Not reached")
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\keyword{survival}

