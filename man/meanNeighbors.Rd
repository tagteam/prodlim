\name{meanNeighbors}
\alias{meanNeighbors}
\title{Helper function to obtain running means for prodlim objects.}
\usage{
  meanNeighbors(x, y, ...)
}
\arguments{
  \item{x}{Object of class \code{"neighborhood"}.}

  \item{y}{Vector of numeric values.}

  \item{\dots}{Not used.}
}
\description{
  Compute average values of a variable according to
  neighborhoods.
}
\examples{
meanNeighbors(x=1:10,y=c(1,10,100,1000,1001,1001,1001,1002,1002,1002))
}
\author{
  Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
  \code{\link{neighborhood}}
}
\keyword{survival}

