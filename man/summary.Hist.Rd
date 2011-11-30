\name{summary.Hist}
\alias{summary.Hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of event histories}
\description{
  Describe events and censoring patterns of an event history. 
}
\usage{
\method{summary}{Hist}(object, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object with class `Hist' derived with \code{\link{Hist}}}
  \item{verbose}{Logical. If FALSE any printing is supressed. }
  \item{\dots}{Not used}
}
% \details{
% }
\value{
  \code{NULL} for survival and competing risk models.
  For other multi-state models, it is
   a list with the following entries:
   \item{states}{the states of the model}
   \item{transitions}{the transitions between the states}
   \item{trans.frame}{a data.frame with the from and to states of the transitions}
}
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk} }
% \note{ ~~further notes~~ 
%  ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{\code{\link{Hist}}, \code{\link{plot.Hist}}}
\examples{
icensFrame <- data.frame(L=c(1,1,3,4,6),R=c(2,NA,3,6,9),event=c(1,1,1,2,2))
with(icensFrame,summary(Hist(time=list(L,R))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
