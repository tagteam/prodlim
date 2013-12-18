\name{atRisk}
\alias{atRisk}
\title{Drawing numbers of subjects at-risk of experiencing an event below
Kaplan-Meier and Aalen-Johannsen plots.}
\usage{
  atRisk(x, newdata, times, line, col, interspace, cex,
    labels, title = "", pos, adj, dist,
    adjust.labels = TRUE, ...)
}
\arguments{
  \item{x}{an object of class `prodlim' as returned by the
  \code{prodlim} function.}

  \item{newdata}{see \code{plot.prodim}}

  \item{times}{Where to compute the atrisk numbers.}

  \item{line}{Distance of the atrisk numbers from the inner
  plot.}

  \item{col}{The color of the text.}

  \item{interspace}{Distance between rows of atrisk
  numbers.}

  \item{cex}{Passed on to \code{mtext} for both atrisk
  numbers and labels.}

  \item{labels}{Labels for the at-risk rows.}

  \item{title}{Title for the at-risk labels}

  \item{pos}{The value is passed on to the \code{mtext}
  argument \code{at} for the labels (not the atriks
  numbers).}

  \item{adj}{Passed on to \code{mtext} for the labels (not
  the atriks numbers).}

  \item{dist}{If \code{line} is missing, the distance of
  the upper most atrisk row from the inner plotting region:
  par()$mgp[2].}

  \item{adjust.labels}{If \code{TRUE} the labels are left
  adjusted.}

  \item{...}{Further arguments that are passed to the
  function \code{mtext}.}
}
\value{
  Nil
}
\description{
  This function is invoked and controlled by
  \code{plot.prodlim}.
}
\details{
  This function should not be called directly. The
  arguments can be specified as \code{atRisk.arg} in the
  call to \code{plot.prodim}.
}
\author{
  Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
  \code{\link{plot.prodlim}}, \code{\link{confInt}},
  \code{\link{markTime}}
}
\keyword{survival}

