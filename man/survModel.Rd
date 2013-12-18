\name{survModel}
\alias{survModel}
\title{Survival model for simulation}
\usage{
  survModel()
}
\value{
  A structural equation model initialized with three
  variables: the latent event time, the latent right
  censored time, and the observed right censored event
  time.
}
\description{
  Create a survival model to simulate a right censored
  event time data without covariates
}
\details{
  This function requires the \code{lava} package.
}
\author{
  Thomas A. Gerds <tag@biostat.ku.dk>
}

