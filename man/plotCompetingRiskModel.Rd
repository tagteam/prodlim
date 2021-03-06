% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCompetingRiskModel.R
\name{plotCompetingRiskModel}
\alias{plotCompetingRiskModel}
\title{Plotting a competing-risk-model.}
\usage{
plotCompetingRiskModel(stateLabels, horizontal = TRUE, ...)
}
\arguments{
\item{stateLabels}{Labels for the boxes.}

\item{horizontal}{The orientation of the plot.}

\item{\dots}{Arguments passed to \code{\link{plot.Hist}}.}
}
\description{
Plotting a competing-risk-model.
}
\examples{

plotCompetingRiskModel()
plotCompetingRiskModel(labels=c("a","b"))
plotCompetingRiskModel(labels=c("a","b","c"))

}
\seealso{
\code{\link{plotIllnessDeathModel}}, \code{\link{plot.Hist}}
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\keyword{survival}
