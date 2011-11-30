\name{NN}
\alias{NN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dummy function}
\description{
  This is a special function used in the context of the conditional
  product limit estimator. It is similar to the function strata and can be used in
  a formula to force a stratified analysis in  
   nearest neighborhoods defined by the values of a continuous predictor.
}
\usage{
 NN(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A continuous predictor}
}
\value{
  The argument is passed on as it is: NN(x)=x
}
\author{Thomas Gerds }
\seealso{\code{\link{prodlim}}}
\examples{
\dontrun{
require(survival)
data(pbc)
prodlim(Hist(time,status)~NN(age),data=pbc)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
