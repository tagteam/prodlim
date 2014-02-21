\name{predict.prodlim}
\alias{predict.prodlim}
\alias{predictCuminc}
\alias{predictSurv}
\title{Predicting event probabilities from product limit estimates}
\usage{
  \method{predict}{prodlim} (object, times, newdata,
    level.chaos = 1, type = c("surv", "cuminc", "list"),
    mode = "list", bytime = FALSE, cause = 1, ...)
}
\arguments{
  \item{object}{A fitted object of class "prodlim".}

  \item{times}{Vector of times at which to return the
  estimated probabilities.}

  \item{newdata}{A data frame with the same variable names
  as those that appear on the right hand side of the
  'prodlim' formula.  If there are covariates this argument
  is required.}

  \item{level.chaos}{Integer specifying the sorting of the
  output: `0' sort by time and newdata; `1' only by time;
  `2' no sorting at all}

  \item{type}{Choice between "surv","cuminc","list":

  "surv": predict survival probabilities only survival
  models

  "cuminc": predict cumulative incidences only competing
  risk models

  "list": find the indices corresponding to times and
  newdata. See value.

  Defaults to "surv" for two-state models and to "cuminc"
  for competing risk models.}

  \item{mode}{Only for \code{type=="surv"} and
  \code{type=="cuminc"}. Can either be "list" or "matrix".
  For "matrix" the predicted probabilities will be returned
  in matrix form.}

  \item{bytime}{Logical. If TRUE and \code{mode=="matrix"}
  the matrix with predicted probabilities will have a
  column for each time and a row for each newdata. Only
  when \code{object$covariate.type>1} and more than one
  time is given.}

  \item{cause}{The cause for predicting the cause-specific
  cumulative incidence function in competing risk models.}

  \item{\dots}{Only for compatibility reasons.}
}
\value{
  \code{type=="surv"} A list or a matrix with survival
  probabilities for all times and all newdata.

  \code{type=="cuminc"} A list or a matrix with cumulative
  incidences for all times and all newdata.

  \code{type=="list"} A list with the following components:

  \item{times}{The argument \code{times} carried forward}

  \item{predictors}{The relevant part of the argument
  \code{newdata}.} \item{indices}{ A list with the
  following components

  \code{time}: Where to find values corresponding to the
  requested times \code{strata}: Where to find values
  corresponding to the values of the variables in newdata.
  Together time and strata show where to find the predicted
  probabilities.  } \item{dimensions}{ a list with the
  following components: \code{time} : The length of
  \code{times} \code{strata} : The number of rows in
  \code{newdata} \code{names.strata} : Labels for the
  covariate values.  }
}
\description{
  Evaluation of estimated survival or event probabilities
  at given times and covariate constellations.
}
\details{
  Predicted (survival) probabilities are returned that can
  be plotted, summarized and used for inverse of
  probability of censoring weighting.
}
\examples{
dat <- SimSurv(400)
fit <- prodlim(Hist(time,status)~1,data=dat)

## predict the survival probs at selected times
predict(fit,times=c(10,100,1000))

## works also outside the usual range of the Kaplan-Meier
predict(fit,times=c(-1,0,10,100,1000,10000))

## newdata is required if there are strata
## or neighborhoods (i.e. overlapping strata)
mfit <- prodlim(Hist(time,status)~X1+X2,data=dat)
predict(mfit,times=c(-1,0,10,100,1000,10000),newdata=dat[18:21,])

## this can be requested in matrix form
predict(mfit,times=c(-1,0,10,100,1000,10000),newdata=dat[18:21,],mode="matrix")

## and even transposed
predict(mfit,times=c(-1,0,10,100,1000,10000),newdata=dat[18:21,],mode="matrix",bytime=TRUE)
}
\author{
  Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
  \code{\link{predictSurvIndividual}}
}
\keyword{survival}

