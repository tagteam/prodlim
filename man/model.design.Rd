\name{model.design}
\alias{model.design}
\title{Extract a design matrix and specials from a model.frame}
\usage{
model.design(data, xlev, dropIntercept = FALSE, maxOrder = 1,
  unspecialsDesign = TRUE, specialsFactor = TRUE, specialsDesign = FALSE,
  stripSpecialNames = TRUE)
}
\arguments{
  \item{data}{A model.frame}

  \item{xlev}{Passed to model.matrix}

  \item{dropIntercept}{If TRUE drop intercept term from the
  design matrix}

  \item{maxOrder}{An error is produced if special variables
  are involved in interaction terms of order higher than
  max.order.}

  \item{unspecialsDesign}{A logical value: if \code{TRUE}
  apply \code{\link{model.matrix}} to unspecial covariates.
  If \code{FALSE} extract unspecial covariates from data.}

  \item{specialsFactor}{A character vector containing
  special variables which should be coerced into a single
  factor. If \code{TRUE} all specials are treated in this
  way.}

  \item{specialsDesign}{A character vector containing
  special variables which should be transformed into a
  design matrix via \code{\link{model.matrix}}.  If
  \code{TRUE} all specials are treated in this way.}

  \item{stripSpecialNames}{If TRUE strip the special from
  variable name, i.e., use X instead of strata(X).}
}
\value{
A list which contains - the design matrix with the levels
of the variables stored in attribute 'levels' - separate
data.frames which contain the values of the special
variables.
}
\description{
Extract design matrix and data specials from a model.frame
}
\details{
The function \code{untangle.specials} of the survival
function does a similar job.
}
\examples{
f <- formula(y~x+ID(z))
set.seed(8)
d <- data.frame(y=rnorm(5),x=factor(c("a","b","b","a","c")),z=c(2,2,7,7,7))
ID <- function(x)x
t <- terms(f,special="ID",data=d)
m <- model.frame(t,d)
md <- model.design(m,specialsFactor=TRUE)
md
md <- model.design(m,specialsFactor=TRUE,unspecialsDesign=FALSE)
md

# special function with argument
treat <- function(x,...) x
f2 <- formula(y~x+treat(z,arg=2)+treat(u,arg=-1))
set.seed(8)
d <- data.frame(y=rnorm(5),u=1:5,x=factor(c("a","b","b","a","c")),z=c(2,2,7,7,7))
t2 <- terms(f2,special="treat",data=d)
m2 <- model.frame(t2,d)
md2 <- model.design(m2,specialsFactor=TRUE)

library(survival)
data(pbc)
tt <- terms(Surv(time,status!=0)~factor(edema)*age+strata(I(log(bili)>1))+strata(sex)+cluster(id),
special=c("strata","cluster"),data=pbc[1:10,])
dd <- model.frame(tt,data=pbc[1:10,])
model.design(dd)
model.design(dd,dropIntercept=TRUE)
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
\code{\link{EventHistory.frame}} model.frame terms
model.matrix .getXlevels
}

