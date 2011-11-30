\name{clusterTestData}
\alias{clusterTestData}
\docType{data}
\title{
  Data set with a clustered data structure,  for testing the product limit method.
}
\description{
  Data set with a clustered data structure,  for testing the product limit method.
}
\usage{data(clusterTestData)}
\format{
  A data frame with 42 observations on the following 4 variables.
  \describe{
    \item{\code{midtimeX}}{event time}
    \item{\code{eventX}}{event status}
    \item{\code{patientid}}{patient number}
    \item{\code{AnyCrownFracture}}{a binary predictor}
  }
}
\examples{
data(clusterTestData)
## maybe str(clusterTestData) ; plot(clusterTestData) ...
}
\keyword{datasets}
