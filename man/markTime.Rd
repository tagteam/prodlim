\name{markTime}
\Rdversion{1.1}
\alias{markTime}
\title{
  Marking product-limit plots at the censored times.
}
\description{
  This function is invoked and controlled by \code{plot.prodlim}.
}
\usage{
markTime(x, times, nlost, pch, col,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The values of the curves at \code{times}.
  }
  \item{times}{
    The times where there curves are plotted.
  }
  \item{nlost}{
    The number of subjects lost to follow-up (censored) at \code{times}.
  }
  \item{pch}{
    The symbol used to mark the curves.
  }
  \item{col}{
    The color of the symbols.
  }
  \item{...}{
    Arguments passed to \code{points}.
  }
}
\details{
    This function should not be called directly. The arguments can be
  specified as \code{atRisk.arg} in the call to \code{plot.prodim}.
}
\value{
  Nil
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}

\seealso{
  \code{\link{plot.prodlim}}, \code{\link{confInt}}, \code{\link{atRisk}}
}
\keyword{survival}