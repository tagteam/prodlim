\name{NN}
\alias{NN}
\title{Dummy function}
\usage{
NN(x)
}
\arguments{
  \item{x}{A continuous predictor variable}
}
\value{
The argument is passed on as it is: NN(x)=x
}
\description{
This is a special function used in the context of the
conditional product limit estimator. It is similar to the
function strata and can be used in a formula to force a
stratified analysis in nearest neighborhoods defined by the
values of a continuous predictor.
}
\examples{
dat <- SimSurv(20)
prodlim(Hist(time,status)~NN(X2),data=dat)
}
\author{
Thomas Gerds
}
\seealso{
\code{\link{prodlim}}
}
\keyword{nonparametric}

