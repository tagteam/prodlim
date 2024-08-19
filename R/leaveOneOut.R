#' Compute leave-one-out estimates
#'
#' This function is the work-horse for \code{jackknife}
#' @title Compute jackknife pseudo values.
#' @aliases leaveOneOut leaveOneOut.survival leaveOneOut.competing.risks 
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{jackknife}}
#' 
#' @param object Object of class \code{"prodlim"}.
#' @param times time points at which to compute leave-one-out
#' event/survival probabilities.
#' @param cause Character (other classes are converted with \code{as.character}).
#' For competing risks the cause of interest. 
#' @param lag For survival models only. If \code{TRUE} lag the result, i.e. compute
#' S(t-) instead of S(t).
#' @param ... not used
#' @export
leaveOneOut <- function(object,times,cause,lag=FALSE,...){
    if (missing(cause)) cause <- attr(object$model.response,which="states")[[1]]
    else
        cause <- checkCauses(cause,object)
    if (object$model=="survival")
        leaveOneOut.survival(object=object,times=times,lag=lag,...)
    else if (object$model=="competing.risks")
        leaveOneOut.competing.risks(object=object,times=times,cause=cause,...)
    else stop("No method for jackknifing this object.")
}

#' @export 
leaveOneOut.survival <- function(object,times,lag=FALSE,...){
    stopifnot(object$covariate.type==1)
    mr <- object$model.response
    time <- object$time
    Y <- object$n.risk
    D <- object$n.event
    Y <- Y[D>0]
    time <- time[D>0]
    D <- D[D>0]
    NU <- length(time)
    obstimes <- mr[,"time"]
    status <- mr[,"status"]
    N <- length(obstimes)
    S <- predict(object,times=time,newdata=mr)
    ## idea: find the at-risk set for pseudo-value k by
    ##       subtracting 1 in the period where subj k is
    ##       at risk. need the position of obstime.k in time ...
    found <- sindex(jump.times=time,eval.times=times)
    loo <- .C("loo_surv",
              Y = as.double(Y),
              D=as.double(D),
              time=as.double(time),
              obsT=as.double(obstimes),
              status=as.double(status),
              S=double(NU),
              loo=double(N*length(times)),
              N=as.integer(N),
              NT=as.integer(NU),
              NP=as.integer(length(times)),
              pos=as.integer(found),
              lag=as.integer(lag),
              PACKAGE="prodlim")$loo
    matrix(loo,ncol=length(times),byrow=FALSE)
}

#' @export 
leaveOneOut.competing.risks <- function(object,times,cause,...){
    stopifnot(object$covariate.type==1)
    mr <- object$model.response
    states <- attr(mr,"states")
    if (missing(cause)) {
        C <- 1
        cause <- states[1]
    }
    else{
        C <- match(cause,states,nomatch=0)
        if (length(C)>1 || C==0) stop("Cause must match exactly one of the names of object$n.event.")
    }
    D <- object$n.event[[C]]
    Dall <- Reduce("+",object$n.event)
    #  it is sufficient to consider time points where events occur
    time <- object$time[Dall>0]
    Y <- object$n.risk[Dall>0]
    D <- D[Dall>0]
    Dall <- Dall[Dall>0]
    NU <- length(time)
    obstimes <- mr[,"time"]
    status <- mr[,"status"]
    E <- getEvent(mr)
    N <- length(obstimes)
    ## idea: see leaveOneOut.survival
    found <- sindex(jump.times=time,eval.times=times)
    ## print(cbind(D,Dall,time,Y))
    loo <- .C("loo_comprisk",
              Y = as.double(Y),
              D=as.double(D),
              Dall=as.double(Dall),
              time=as.double(time),
              obsT=as.double(obstimes),
              status=as.double(status),
              event=as.double(status*(E==cause)),
              S=double(NU),
              F=double(NU),
              loo=double(N*length(times)),
              N=as.integer(N),
              NT=as.integer(NU),
              NP=as.integer(length(times)),
              pos=as.integer(found),
              PACKAGE="prodlim")$loo
    matrix(loo,ncol=length(times),byrow=FALSE)
}
