\name{summary.prodlim}
\alias{summary.prodlim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary method for prodlim objects.
}
\description{
  Summarizing the result of the product limit method in life-table
  format. Calculates the number of subjects at
  risk and counts events and censored observations at specified times or
  in specified time intervals.
}
\usage{
\method{summary}{prodlim}(object,
		times,
		newdata,
		max.tables = 20,
                surv=TRUE,
                cause,
		intervals = FALSE,
		percent = FALSE,
                showTime=TRUE,
		...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    An object with class `prodlim' derived with \code{\link{prodlim}}}
  
  \item{times}{
    Vector of times at which to return the estimated probabilities.
  }
  
  \item{newdata}{
    A data frame with the same variable names as those that
    appear on the right hand side of the 'prodlim' formula.
    Defaults to \code{object$model.matrix}.
  }
  \item{max.tables}{
    Integer. If \code{newdata} is not given the value of
    \code{max.tables} decides about the maximal number
    of tables to be shown. Defaults to 20.
  }
  \item{surv}{
    Logical. If FALSE report event probabilities 
    instead of survival probabilities. Only available for
    \code{object$model=="survival"}.
  }
\item{cause}{
  The cause for predicting the cause-specific cumulative incidence
  function in competing risk models.
}
  
  \item{intervals}{
    Logical. If TRUE count events and censored in intervals between the values of
    \code{times}.
  }
  \item{percent}{
    Logical. If TRUE all estimated values are multiplied by 100 and thus
    interpretable on a percent scale.
  }
  \item{showTime}{
    Passed to \code{lifetable}
  }
  \item{\dots}{
    Further arguments that are passed to the print function.
  }
}
\details{
  For cluster-correlated data the number of clusters at-risk are 
  are also given. Confidence intervals are displayed when they
  are part of the fitted object. 
}
\value{
  A data.frame with the relevant information.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk} }
% \note{
%   ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{prodlim}}, \code{\link{summary.Hist}}}
\examples{
SurvFrame <- data.frame(time=1:10,status=rbinom(10,1,.5))
 x <- prodlim(formula=Hist(time=time,status!=0)~1,data=SurvFrame)
 summary(x)
 summary(x,times=c(1,5,10,12),percent=TRUE,intervals=TRUE,digits=3)
\dontrun{
library(survival) 
data(pbc)
f <- prodlim(Hist(time,status!=1)~sex,data=pbc)
summary(f)
summary(f,newdata=data.frame(sex=c("m","f")))

fa <- prodlim(Hist(time,status!=1)~sex+age,data=pbc)
summary(fa)

x <- summary(fa,times=1000,newdata=expand.grid(age=c(60,40,50),sex=c("m","f")))
cbind(names(x$table),do.call("rbind",lapply(x$table,round,2)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
