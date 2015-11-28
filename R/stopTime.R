### stopTime.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Nov 28 2015 (10:07) 
## Version: 
## last-updated: Nov 28 2015 (10:29) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' All event times are stopped at a given time point and
##' corresponding events are censored
##'
##' @title Stop the time of an event history object
##' @param object Event history object as obtained with \code{Hist}
##' @param time Time point at which to stop the event history object
##' @return Stopped object
##' @seealso Hist 
##' @examples
##'
##' d <- SimSurv(10)
##' h <- with(d,Hist(time,status))
##' stopTime(h,5)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
stopTime <- function(object,time){
    if (missing(time)) stop("Argument time missing. Need a time point at which to stop the event history.")
    if (length(time)>1) {
        warning("Argument time is a vector. Proceed with the first element.")
        time <- time[[1]]
    }
    stopifnot(class(object)[[1]]=="Hist")
    model <- attr(object,"model")
    if(!(model %in% c("survival","competing.risks")))
        stop(paste("Don't know yet how to stop time of this type of model:",model))
    censored <- object[,"time"] >=time
    object[,"status"][censored] <- 0
    if(model=="competing.risks")
        object[,"event"][censored] <- attr(object,"cens.code")
    object[,"time"][censored] <- time
    attr(object,"stop.time") <- time
    object
}
#----------------------------------------------------------------------
### stopTime.R ends here
