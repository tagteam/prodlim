\name{PercentAxis}
\alias{PercentAxis}
\title{Percentage-labeled axis.}
\usage{
PercentAxis(x, at, ...)
}
\arguments{
  \item{x}{Side of the axis}

  \item{at}{Positions (decimals) at which to label the
  axis.}

  \item{\dots}{Given to \code{axis}.}
}
\description{
Use percentages instead of decimals to label the an axis
with a probability scale .
}
\examples{
plot(0,0,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
  PercentAxis(1,at=seq(0,1,.25))
  PercentAxis(2,at=seq(0,1,.25))
}
\author{
Thomas Alexander Gerds
}
\seealso{
\code{\link{plot.prodlim}}
}
\keyword{survival}

