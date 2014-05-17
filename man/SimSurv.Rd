\name{SimSurv}
\alias{SimSurv}
\title{Simulate survival data}
\usage{
SimSurv(N, ...)
}
\arguments{
  \item{N}{sample size}

  \item{...}{do nothing}
}
\value{
data.frame with simulated data
}
\description{
Simulate right censored survival data with two covariates
X1 and X2, both have effect exp(1) on the hazard of the
unobserved event time.
}
\details{
This function calls \code{survModel}, then adds covariates
and finally calls \code{sim.lvm}.
}
\examples{
SimSurv(10)
}
\author{
Thomas Alexander Gerds
}
\references{
Bender, Augustin & Blettner. Generating survival times to
simulate Cox proportional hazards models. Statistics in
Medicine, 24: 1713-1723, 2005.
}

