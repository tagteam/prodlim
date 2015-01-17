\name{plot.prodlim}
\alias{lines.prodlim}
\alias{plot.prodlim}
\title{Plotting event probabilities over time}
\usage{
\method{plot}{prodlim}(x, type, cause = 1, select, newdata, add = FALSE,
  col, lty, lwd, ylim, xlim, xlab = "Time", ylab, legend = TRUE,
  logrank = FALSE, marktime = FALSE, confint = TRUE, automar,
  atrisk = ifelse(add, FALSE, TRUE), timeOrigin = 0, axes = TRUE,
  background = TRUE, percent = TRUE, minAtrisk = 0, limit = 10, ...)
}
\arguments{
  \item{x}{an object of class `prodlim' as returned by the
  \code{prodlim} function.}

  \item{type}{controls what part of the object is plotted.
  Defaults to \code{"survival"} for the Kaplan-Meier
  estimate of the survival function in two state models and
  to \code{"incidence"} for the Aalen-Johansen estimate of
  the cumulative incidence function in competing risk
  models}

  \item{cause}{determines the cause of the cumulative
  incidence function.  Currently one cause is allowed at a
  time, but you may call the function again with add=TRUE
  to add the lines of the other causes.}

  \item{select}{Select which lines to plot. This can be
  used when there are many strata coming out by default to
  quickly select a subset of the lines (however, perhaps
  better to use newdata argument). Useful also for stacked
  plot when there are many events and one wants to stack
  only a subset of the cumulative incidence functions.}

  \item{newdata}{a data frame containing strata for which
  plotted curves are desired. When omitted element \code{X}
  of object \code{x} is used.}

  \item{add}{if 'TRUE' curves are added to an existing
  plot.}

  \item{col}{color for curves defaults to 1:number(curves)}

  \item{lty}{line type for curves defaults to 1}

  \item{lwd}{line width for all curves}

  \item{ylim}{limits of the y-axis}

  \item{xlim}{limits of the x-axis}

  \item{xlab}{label for the x-axis}

  \item{ylab}{label for the y-axis}

  \item{legend}{if TRUE a legend is plotted by calling the
  function legend.  Optional arguments of the function
  \code{legend} can be given in the form
  \code{legend.x=val} where x is the name of the argument
  and val the desired value. See also Details.}

  \item{logrank}{If TRUE, the logrank p-value will be
  extracted from a call to \code{survdiff} and added to the
  legend. This works only for survival models, i.e.
  Kaplan-Meier with discrete predictors.}

  \item{marktime}{if TRUE the curves are tick-marked at
  right censoring times by invoking the function
  \code{markTime}. Optional arguments of the function
  \code{markTime} can be given in the form
  \code{confint.x=val} as with legend. See also Details.}

  \item{confint}{if TRUE pointwise confidence intervals are
  plotted by invoking the function \code{confInt}. Optional
  arguments of the function \code{confInt} can be given in
  the form \code{confint.x=val} as with legend.  See also
  Details.}

  \item{automar}{If TRUE the function trys to find suitable
  values for the figure margins around the main plotting
  region.}

  \item{atrisk}{if TRUE display numbers of subjects at risk
  by invoking the function \code{atRisk}. Optional
  arguments of the function \code{atRisk} can be given in
  the form \code{atrisk.x=val} as with legend. See also
  Details.}

  \item{timeOrigin}{Start of the time axis}

  \item{axes}{If true axes are drawn. See details.}

  \item{background}{If \code{TRUE} the background color and
  grid color can be controlled using smart arguments
  SmartControl, such as background.bg="yellow" or
  background.bg=c("gray66","gray88").  The following
  defaults are passed to \code{background} by
  \code{plot.prodlim}: horizontal=seq(0,1,.25),
  vertical=NULL, bg="gray77", fg="white".  See
  \code{background} for all arguments, and the examples
  below.}

  \item{percent}{If true the y-axis is labeled in percent.}

  \item{minAtrisk}{Integer. Show the curve only until the
  number at-risk is at least \code{minAtrisk}}

  \item{limit}{When newdata is not specified and the number
  of lines in element \code{X} of object \code{x} exceeds
  limits, only the results for covariate constellations of
  the first, the middle and the last row in \code{X} are
  shown. Otherwise all lines of \code{X} are shown.}

  \item{...}{Parameters that are filtered by
  \code{\link{SmartControl}} and then passed to the
  functions \code{\link{plot}}, \code{\link{legend}},
  \code{\link{axis}}, \code{\link{atRisk}},
  \code{\link{confInt}}, \code{\link{markTime}},
  \code{\link{backGround}}}
}
\value{
The (invisible) object.
}
\description{
Function to plot survival and cumulative incidence curves
against time.
}
\details{
From version 1.1.3 on the arguments legend.args,
atrisk.args, confint.args are obsolete and only available
for backward compatibility. Instead arguments for the
invoked functions \code{atRisk}, \code{legend},
\code{confInt}, \code{markTime}, \code{axis} are simply
specified as \code{atrisk.cex=2}. The specification is not
case sensitive, thus \code{atRisk.cex=2} or
\code{atRISK.cex=2} will have the same effect.  The
function \code{axis} is called twice, and arguments of the
form \code{axis1.labels}, \code{axis1.at} are used for the
time axis whereas \code{axis2.pos}, \code{axis1.labels},
etc. are used for the y-axis.

These arguments are processed via \code{\dots{}} of
\code{plot.prodlim} and inside by using the function
\code{SmartControl}.  Documentation of these arguments can
be found in the help pages of the corresponding functions.
}
\note{
Similar functionality is provided by the function
\code{\link{plot.survfit}} of the survival library
}
\examples{
## simulate right censored data from a two state model
set.seed(100)
dat <- SimSurv(100)
# with(dat,plot(Hist(time,status)))

### marginal Kaplan-Meier estimator
kmfit <- prodlim(Hist(time, status) ~ 1, data = dat)
plot(kmfit)

# change time range
plot(kmfit,xlim=c(0,4))

# change scale of y-axis
plot(kmfit,percent=FALSE)

# change axis label and position of ticks
plot(kmfit,
     xlim=c(0,10),
     axis1.at=seq(0,10,1),
     axis1.labels=0:10,
     xlab="Years",
     axis2.las=2,
     atrisk.at=seq(0,10,2.5),
     atrisk.title="")

# change background color
plot(kmfit,
     xlim=c(0,10),
     confint.citype="shadow",
     col=1,
     axis1.at=0:10,
     axis1.labels=0:10,
     xlab="Years",
     axis2.las=2,
     atrisk.at=seq(0,10,2.5),
     atrisk.title="",
     background=TRUE,
     background.fg="white",
     background.horizontal=seq(0,1,.25/2),
     background.vertical=seq(0,10,2.5),
     background.bg=c("gray88"))

# change type of confidence limits
plot(kmfit,
     xlim=c(0,10),
     confint.citype="dots",
     col=4,
     background=TRUE,
     background.bg=c("white","gray88"),
     background.fg="gray77",
     background.horizontal=seq(0,1,.25/2),
     background.vertical=seq(0,10,2))


### Kaplan-Meier in discrete strata
kmfitX <- prodlim(Hist(time, status) ~ X1, data = dat)
plot(kmfitX)
# move legend
plot(kmfitX,legend.x="bottomleft",atRisk.cex=1.3,
     atrisk.title="No. subjects")

## Control the order of strata
## since version 1.5.1 prodlim does  obey the order of
## factor levels
dat$group <- factor(cut(dat$X2,c(-Inf,0,0.5,Inf)),
                    labels=c("High","Intermediate","Low"))
kmfitG <- prodlim(Hist(time, status) ~ group, data = dat)
plot(kmfitG)

## relevel
dat$group2 <- factor(cut(dat$X2,c(-Inf,0,0.5,Inf)),
                     levels=c("(0.5, Inf]","(0,0.5]","(-Inf,0]"),
                     labels=c("Low","Intermediate","High"))
kmfitG2 <- prodlim(Hist(time, status) ~ group2, data = dat)
plot(kmfitG2)

# add log-rank test to legend
plot(kmfitX,
     atRisk.cex=1.3,
     logrank=TRUE,
     legend.x="topright",
     atrisk.title="at-risk")

# change atrisk labels
plot(kmfitX,
     legend.x="bottomleft",
     atrisk.title="Patients",
     atrisk.cex=0.9,
     atrisk.labels=c("X1=0","X1=1"))

# multiple categorical factors

kmfitXG <- prodlim(Hist(time,status)~X1+group2,data=dat)
plot(kmfitXG,select=1:2)

### Kaplan-Meier in continuous strata
kmfitX2 <- prodlim(Hist(time, status) ~ X2, data = dat)
plot(kmfitX2,xlim=c(0,10))

# specify values of X2 for which to show the curves
plot(kmfitX2,xlim=c(0,10),newdata=data.frame(X2=c(-1.8,0,1.2)))

### Cluster-correlated data
library(survival)
cdat <- cbind(SimSurv(20),patnr=sample(1:5,size=20,replace=TRUE))
kmfitC <- prodlim(Hist(time, status) ~ cluster(patnr), data = cdat)
plot(kmfitC)
plot(kmfitC,atrisk.labels=c("Units","Patients"))

kmfitC2 <- prodlim(Hist(time, status) ~ X1+cluster(patnr), data = cdat)
plot(kmfitC2)
plot(kmfitC2,atrisk.labels=c("Teeth","Patients","Teeth","Patients"),
     atrisk.col=c(1,1,2,2))

## simulate right censored data from a competing risk model
datCR <- SimCompRisk(100)
with(datCR,plot(Hist(time,event)))

### marginal Aalen-Johansen estimator
ajfit <- prodlim(Hist(time, event) ~ 1, data = datCR)
plot(ajfit) # same as plot(ajfit,cause=1)

# cause 2
plot(ajfit,cause=2)

# both in one
plot(ajfit,cause=1)
plot(ajfit,cause=2,add=TRUE,col=2)

### stacked plot

plot(ajfit,cause="stacked",select=2)

### conditional Aalen-Johansen estimator
ajfitX1 <- prodlim(Hist(time, event) ~ X1, data = datCR)
plot(ajfitX1)

ajfitX <- prodlim(Hist(time, event) ~ X1+X2, data = datCR)
plot(ajfitX,newdata=data.frame(X1=c(1,1,0),X2=c(4,10,10)))
plot(ajfitX,newdata=data.frame(X1=c(1,1,0),X2=c(4,10,10)),cause=2)

## stacked plot

plot(ajfitX,
     newdata=data.frame(X1=0,X2=0.1),
     cause="stacked",
     legend.title="X1=0,X2=0.1",
     legend.legend=paste("cause:",getStates(ajfitX$model.response)),
     plot.main="Subject specific stacked plot")
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
\code{\link{plot}}, \code{\link{legend}},
\code{\link{axis}},
\code{\link{prodlim}},\code{\link{plot.Hist}},\code{\link{summary.prodlim}},
\code{\link{neighborhood}}, \code{\link{atRisk}},
\code{\link{confInt}}, \code{\link{markTime}},
\code{\link{backGround}}
}
\keyword{survival}

