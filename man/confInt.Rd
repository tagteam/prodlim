\name{confInt}
\alias{confInt}
\title{Add point-wise confidence limits to the graphs of Kaplan-Meier and
Aalen-Johannsen estimates of survival and cumulative incidence.}
\usage{
confInt(x, times, newdata, type, citype, cause, col, lty, lwd, density = 55,
  ...)
}
\arguments{
  \item{x}{an object of class `prodlim' as returned by the
  \code{prodlim} function.}

  \item{times}{where to compute point-wise confidence
  limits}

  \item{newdata}{see \code{plot.prodim}}

  \item{type}{Either \code{"cuminc"} or \code{"survival"}
  passed to summary.prodlim as
  \code{surv=ifelse(type=="cuminc",FALSE,TRUE)}.}

  \item{citype}{If \code{"shadow"} then confidence limits
  are drawn as colored shadows.  Otherwise, dotted lines
  are used to show the upper and lower confidence limits.}

  \item{cause}{see \code{plot.prodim}}

  \item{col}{the colour of the lines.}

  \item{lty}{the line type of the lines.}

  \item{lwd}{the line thickness of the lines.}

  \item{density}{For \code{citype="shadow"}, the density of
  the shade. Default is 55 percent.}

  \item{\dots}{Further arguments that are passed to the
  function \code{segments} if \code{type=="bars"} and to
  \code{lines} else.}
}
\value{
Nil
}
\description{
This function is invoked and controlled by
\code{plot.prodlim}.
}
\details{
This function should not be called directly. The arguments
can be specified as \code{Confint.arg} in the call to
\code{plot.prodim}.
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
\code{\link{plot.prodlim}}, \code{\link{atRisk}},
\code{\link{markTime}}
}
\keyword{survival}

