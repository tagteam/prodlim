\name{parseSpecialNames}
\alias{parseSpecialNames}
\title{Parse special terms}
\usage{
parseSpecialNames(x, specials, specialArgumentNames)
}
\arguments{
  \item{x}{Vector of character strings, usually the column
  names of the design matrix obtained with
  \code{\link{model.design}}.}

  \item{specials}{A vector with character strings providing
  the names of the special arguments.}

  \item{specialArgumentNames}{A named list with one element
  for each special of argument}
}
\value{
A named list of parsed arguments. The names of the list are
the variable names.
}
\description{
Extract from a vector of character strings the names of
special functions and auxiliary arguments
}
\details{
The names of the special functions may not be nested, i.e.,
c(treat, treatment) will not work.
}
\examples{
parseSpecialNames("treat(Z)",specials="treat")
parseSpecialNames("treat(Z,u=2)",specials="treat",specialArgumentNames=list("treat"="u"))
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
model.design
}

