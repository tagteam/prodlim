\name{predictSurvIndividual}
\alias{predictSurvIndividual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict individual survival probabilities}
\description{
  Function to extract the predicted probabilities at the
  individual event times that have been used for fitting a prodlim
  object.
}
\usage{
predictSurvIndividual(object, lag = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A fitted object of class "prodlim".
  }
  \item{lag}{
    Integer. `0' means predictions at the individual
    times, 1 means just before the individual times, etc.
  }
}
% \details{
% }
\value{
  A vector of survival probabilities.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk} }
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{predict.prodlim}},\code{\link{predictSurv}},}
\examples{
 SurvFrame <- data.frame(time=1:10,status=rbinom(10,1,.5))
 x <- prodlim(formula=Hist(time=time,status!=0)~1,data=SurvFrame)
 predictSurvIndividual(x,lag=1)
}
\keyword{ survival }% at least one, from doc/KEYWORDS
