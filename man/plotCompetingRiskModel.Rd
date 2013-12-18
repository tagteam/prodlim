\name{plotCompetingRiskModel}
\alias{plotCompetingRiskModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plotting a competing-risk-model.
}
\description{
  Plotting a competing-risk-model.
}
\usage{
plotCompetingRiskModel(stateLabels, horizontal = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stateLabels}{
    Labels for the boxes.
  }
  \item{horizontal}{
    The orientation of the plot.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.Hist}}.
 }
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
\code{\link{plotIllnessDeathModel}}, \code{\link{plot.Hist}}
}
\examples{
plotCompetingRiskModel()
plotCompetingRiskModel(labels=c("a","b"))
plotCompetingRiskModel(labels=c("a","b","c"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
