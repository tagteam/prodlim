\name{PercentAxis}
\alias{PercentAxis}
\title{
  Percentage-labeled axis. 
}
\description{
  Use percentages instead of decimals to label the an axis with a
  probability scale .
}
\usage{
PercentAxis(x, at, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Side of the axis
}
  \item{at}{
    Positions (decimals) at which to label the axis.
  }
  \item{\dots}{
    Given to \code{axis}.
  }
}
\author{
  Thomas Alexander Gerds 
}
\seealso{
  \code{\link{plot.prodlim}}
}
\examples{
  plot(0,0,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
  PercentAxis(1,at=seq(0,1,.25))
  PercentAxis(2,at=seq(0,1,.25))
}
\keyword{ survival }
