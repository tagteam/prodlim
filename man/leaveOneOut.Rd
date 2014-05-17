\name{leaveOneOut}
\alias{leaveOneOut}
\alias{leaveOneOut.competing.risks}
\alias{leaveOneOut.survival}
\title{Compute jackknife pseudo values.}
\usage{
leaveOneOut(object, times, cause, lag = FALSE, ...)
}
\arguments{
  \item{object}{Object of class \code{"prodlim"}.}

  \item{times}{time points at which to compute
  leave-one-out event/survival probabilities.}

  \item{cause}{For competing risks the cause of interest.}

  \item{lag}{For survival models only. If \code{TRUE} lag
  the result, i.e. compute S(t-) instead of S(t).}

  \item{...}{not used}
}
\description{
Compute leave-one-out estimates
}
\details{
This function is the work-horse for \code{jackknife}
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
\code{\link{jackknife}}
}

