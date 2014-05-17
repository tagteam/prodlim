\name{sindex}
\alias{sindex}
\title{Index for evaluation of step functions.}
\usage{
sindex(jump.times, eval.times, comp = "smaller", strict = FALSE)
}
\arguments{
  \item{jump.times}{Numeric vector: e.g.\ the unique jump
  times of a step function.}

  \item{eval.times}{Numeric vector: e.g.\ the times where
  the step function should be evaluated}

  \item{strict}{If TRUE make the comparison of jump times
  and eval times strict}

  \item{comp}{If "greater" count the number of jump times
  that are greater (greater or equal when strict==FALSE)
  than the eval times}
}
\value{
Index of the same length as \code{eval.times} containing
the numbers of the \code{jump.times} that are smaller than
or equal to \code{eval.times}.
}
\description{
Returns an index of positions. Intended for evaluating a
step function at selected times. The function counts how
many elements of a vector, e.g. the jump times of the step
function, are smaller or equal to the elements in a second
vector, e.g. the times where the step function should be
evaluated.
}
\details{
If all \code{jump.times} are greater than a particular
\code{eval.time} the sindex returns \code{0}. This must be
considered when sindex is used for subsetting, see the
Examples below.
}
\examples{
test <- list(time = c(1, 1,5,5,2,7,9),
		status = c(1,0,1,0,1,1,0))
fit <- prodlim(Hist(time,status)~1,data=test)
jtimes <- fit$time
etimes <- c(0,.5,2,8,10)
fit$surv
c(1,fit$surv)[1+sindex(jtimes,etimes)]
}
\author{
Thomas A. Gerds \email{tag@biostat.ku.dk}
}
\keyword{misc}

