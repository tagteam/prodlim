\name{sindex}
\alias{sindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Index for evaluation of step functions.
}
\description{Returns an index of positions. Intended for evaluating a step
  function at selected times. The function counts how many elements of a
  vector, e.g. the jump times of the step function, 
  are smaller or equal to the elements in a second vector, e.g. the times
  where the step function should be evaluated.  
}
\usage{
sindex(jump.times, eval.times,comp="smaller",strict=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jump.times}{Numeric vector: e.g.\ the unique jump times of a step function.}
  \item{eval.times}{Numeric vector: e.g.\ the times where the step
    function should be evaluated}
  \item{strict}{If TRUE make the comparison of jump times and eval times
    strict}
  \item{comp}{If "greater" count the number of jump times that are
    greater (greater or equal when strict==FALSE) than the eval times
  }
}
\details{
  If all \code{jump.times} are greater than a particular
  \code{eval.time} the sindex returns \code{0}. This must
  be considered when sindex is used for subsetting, see the
  Examples below. 
}
\value{
  Index of the same length as \code{eval.times} containing the numbers
  of the \code{jump.times} that are smaller than or equal to 
  \code{eval.times}.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
% \note{ ~~further notes~~ }
%  ~Make other sections like Warning with \section{Warning }{....} ~
% \seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{

test <- list(time = c(1, 1,5,5,2,7,9),
		status = c(1,0,1,0,1,1,0))
fit <- prodlim(Hist(time,status)~1,data=test)
jtimes <- fit$time
etimes <- c(0,.5,2,8,10)
fit$surv
c(1,fit$surv)[1+sindex(jtimes,etimes)]

}
\keyword{misc}% at least one, from doc/KEYWORDS
