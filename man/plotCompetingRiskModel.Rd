\name{plotCompetingRiskModel}
\alias{plotCompetingRiskModel}
\title{Plotting a competing-risk-model.}
\usage{
plotCompetingRiskModel(stateLabels, horizontal = TRUE, ...)
}
\arguments{
  \item{stateLabels}{Labels for the boxes.}

  \item{horizontal}{The orientation of the plot.}

  \item{\dots}{Arguments passed to
  \code{\link{plot.Hist}}.}
}
\description{
Plotting a competing-risk-model.
}
\examples{
plotCompetingRiskModel()
plotCompetingRiskModel(labels=c("a","b"))
plotCompetingRiskModel(labels=c("a","b","c"))
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
\code{\link{plotIllnessDeathModel}},
\code{\link{plot.Hist}}
}
\keyword{survival}

