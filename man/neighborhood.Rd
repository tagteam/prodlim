\name{neighborhood}
\alias{neighborhood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Symmetric neighborhoods for kernel smoothing}
\description{
  Nearest neighborhoods for the values of a continuous predictor. The
  result is used for the conditional Kaplan-Meier estimator
  and other conditional product limit estimators.
}
\usage{
neighborhood(x, bandwidth = NULL, kernel = "box")
}
\arguments{
  \item{x}{Numeric vector -- typically the observations of a continuous
    random variate.}
  \item{bandwidth}{Controls the distance between neighbors in a
    neighborhood. It can be a decimal, i.e.\ the bandwidth, or the string
    `"smooth"', in which case \code{N^{-1/4}} is used, \code{N} being the
    sample size,
    or \code{NULL} in which case the \code{\link{dpik}}
    function of the package KernSmooth is used to find the optimal
    bandwidth.}
  \item{kernel}{Only the rectangular kernel ("box") is implemented.}
}
% \details{
% }
\value{
  An object of class 'neighborhood'.
  The value is a list that includes the unique values of `x'
  (\code{values}) for which a neighborhood, consisting of the nearest neighbors, is defined by
  the first neighbor (\code{first.nbh}) of the usually very long
  vector \code{neighbors}
  and the size of the neighborhood (\code{size.nbh}). 
  
  Further values are the arguments \code{bandwidth}, \code{kernel}, the total sample
  size \code{n} and the number of unique values \code{nu}.}

\references{ 
  Stute, W. "Asymptotic Normality of Nearest Neighbor Regression Function
  Estimates", \emph{The Annals of Statistics},
  1984,12,917--926.
}

\author{Thomas Gerds}
% \note{ ~~further notes~~ }

% ~Make other sections like Warning with
% \section{Warning }{....} ~

\seealso{\code{\link{dpik}}, \code{\link{prodlim}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\dontrun{
library(survival)
data(pbc)
neighborhood(pbc$age)
}}
\keyword{smooth}% at least one, from doc/KEYWORDS
