\name{getStates}
\alias{getStates}
\title{States of a multi-state model}
\usage{
  getStates(object, ...)
}
\arguments{
  \item{object}{}

  \item{...}{}
}
\value{
  A character vector with the states of the model.
}
\description{
  Extract the states of a multi-state model
}
\details{
  Applying this function to the fit of prodlim means to
  apply it to \code{fit$model.response}.
}
\author{
  Thomas A. Gerds
}

