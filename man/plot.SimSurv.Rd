\name{plot.SimSurv}
\alias{plot.SimSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot simulated survival data }
\description{
  Survival data are simulated with \code{SimSurv}. 
  Then this function applies the LinksKaplan-Meier estimator 
  -- possibly in strata defined by covariates -- and plots the result.
}
\usage{
\method{plot}{SimSurv}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data.frame as obtained with \code{SimSurv} }
  \item{\dots}{passed to \code{plot.prodlim}}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
% \note{ ~~further notes~~ 
%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{SimSurv}}, \code{\link{plot.prodlim}}}
\examples{
sdat=SimSurv(100)
plot(sdat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
