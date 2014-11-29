# {{{ header
#' Summary method for prodlim objects.
#' 
#' Summarizing the result of the product limit method in life-table format.
#' Calculates the number of subjects at risk and counts events and censored
#' observations at specified times or in specified time intervals.
#' 
#' For cluster-correlated data the number of clusters at-risk are are also
#' given. Confidence intervals are displayed when they are part of the fitted
#' object.
#' 
#' @param object
#' 
#' An object with class `prodlim' derived with \code{\link{prodlim}}
#' @param times Vector of times at which to return the estimated probabilities.
#' @param newdata A data frame with the same variable names as those that
#' appear on the right hand side of the 'prodlim' formula.  Defaults to
#' \code{object$X}.
#' @param max.tables Integer. If \code{newdata} is not given the value of
#' \code{max.tables} decides about the maximal number of tables to be shown.
#' Defaults to 20.
#' @param surv Logical. If FALSE report event probabilities instead of survival
#' probabilities. Only available for \code{object$model=="survival"}.
#' @param cause The cause for predicting the cause-specific cumulative
#' incidence function in competing risk models.
#' @param intervals Logical. If TRUE count events and censored in intervals
#' between the values of \code{times}.
#' @param percent Logical. If TRUE all estimated values are multiplied by 100
#' and thus interpretable on a percent scale.
#' @param showTime If \code{TRUE} evaluation times are put into a column of the
#' output table, otherwise evaluation times are shown as rownames.
#' @param \dots Further arguments that are passed to the print function.
#' @return A data.frame with the relevant information.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{prodlim}}, \code{\link{summary.Hist}}
#' @keywords survival
#' @examples
#' 
#' SurvFrame <- data.frame(time=1:10,status=rbinom(10,1,.5))
#'  x <- prodlim(formula=Hist(time=time,status!=0)~1,data=SurvFrame)
#'  summary(x)
#'  summary(x,times=c(1,5,10,12),percent=TRUE,intervals=TRUE,digits=3)
#' \donttest{
#' library(survival) 
#' data(pbc)
#' f <- prodlim(Hist(time,status!=1)~sex,data=pbc)
#' summary(f)
#' summary(f,newdata=data.frame(sex=c("m","f")))
#' 
#' summary(f,newdata=data.frame(sex=c("m","f","f","m")))
#' 
#' fa <- prodlim(Hist(time,status!=1)~sex+age,data=pbc)
#' summary(fa)
#' 
#' x <- summary(fa,times=1000,newdata=expand.grid(age=c(60,40,50),sex=c("m","f")))
#' cbind(names(x$table),do.call("rbind",lapply(x$table,round,2)))
#' 
#' }
#'
#' @S3method summary prodlim
#' @method summary prodlim
summary.prodlim <- function(object,
                            times,
                            newdata,
                            max.tables=20,
                            surv=TRUE,
                            cause,
                            intervals=FALSE,
                            percent=FALSE,
                            showTime=TRUE,
                            ...) {
    # }}}
    # {{{  classify the situation
    cens.type <- object$cens.type         # uncensored, right or interval censored
    model <- object$model                 # survival, competing risks or multi-state
    ## cluster <- object$clustervar          # clustered data?
    cotype <- object$covariate.type       # no, discrete, continuous or both
    # }}}
    # {{{  times
    jump.times <- object$time
    if (missing(times) && (length(times <- jump.times) > 50)) 
        times <- quantile(sort(unique(jump.times)))
    times <- sort(unique(times))
    if (any(times>max(jump.times)))
        warning(call.=TRUE,
                immediate.=TRUE,
                paste("\n","Time(s) ",paste(times[times>max(jump.times)],collapse=", "),
                      " are beyond the maximal follow-up time ",max(jump.times),"\n"))
    ntimes <- length(times)
    # }}}
    # {{{ interval-censored
    if (cens.type=="intervalCensored"){
        ltab <- data.frame(time=paste("(",paste(signif(object$time[1,],2),
                               signif(object$time[2,],2),
                               sep="-"),"]",sep=""),
                           n.risk=signif(object$n.risk,2),
                           n.event=signif(object$n.event,2),
                           ##    n.lost=object$n.lost,
                           surv=object$surv)
    }
    else{
        # }}}
        # {{{ with covariates
        if (cotype>1){
            if (missing(newdata) || length(newdata)==0){
                X <- object$X
                if (NROW(X)>max.tables){
                    warning(call.=TRUE,immediate.=TRUE,paste("\nLife tables are available for",
                                           NROW(X),
                                           "different covariate constellations.\n",
                                           "Shown are the table corresponding to the first row in object$X,",
                                           "corresponding to the middle row (median of the number of rows in object$X) ",
                                           "and corresponding to the last row in object$X ...\n",
                                           "to see more tables use arguments `newdata' and `max.tables'\n"))
                    X <- X[c(1,round(median(1:NROW(X))),NROW(X)),,drop=FALSE]
                }
            } else{
                X <- unique.data.frame(newdata)
                if (NROW(X) < NROW(newdata))
                    warning("Returned is only one summary for each unique value in newdata.")
            }
        } else {
            X <- NULL
        }
        if (model=="survival") {
            stats <- list(c("surv",1),c("se.surv",0))
            if (!is.null(object$conf.int))
                stats <- c(stats,list(c("lower",0),c("upper",1)))
            if (surv==FALSE){
                object$cuminc <- 1-object$surv
                object$se.cuminc <- object$se.surv
                cuminc.upper <- 1-object$lower
                cuminc.lower <- 1-object$upper
                object$lower <- cuminc.lower
                object$upper <- cuminc.upper
                stats <- list(c("cuminc",0),c("se.cuminc",0))
                if (!is.null(object$conf.int))
                    stats <- c(stats,list(c("lower",0),c("upper",1)))
            }
        }
        if (model=="competing.risks"){
            stats <- list(c("cuminc",0),c("se.cuminc",0))
            if (!is.null(object$conf.int))
                stats <- c(stats,list(c("lower",0),c("upper",0)))
        }
        ltab <- lifeTab(object=object,
                        times=times,
                        newdata=X,
                        stats=stats,
                        intervals=intervals,
                        percent=percent,
                        showTime=showTime)
        if (model=="competing.risks"){
            if (!missing(cause)){
                if (is.numeric(cause) && !is.numeric(names(ltab)))
                    cause <- attr(object$model.response,"states")[cause]
            } else{ ## show all causes
                cause <- attr(object$model.response,"states")
            }
            # found <- match(cause,attr(object$model.response,"states"),nomatch=FALSE))
            Found <- match(cause,names(ltab))
            if (all(Found)) ltab <- ltab[Found]
            else warning("could not find requested causes in attributes of object$mode.response")
        }
    }
    # }}}
    # {{{ output
    if (is.list(ltab)) {
        tab <- lapply(ltab,function(x){
            if (is.list(x)) {
                lapply(x,data.frame) }
            else{
                data.frame(as.matrix(x))
            }
        })
    }
    else{
        tab <- data.frame(ltab)
    }
    if (model=="competing.risks")
        out <- list(table=ltab,cause=cause)
    else
        out <- list(table=ltab)
    out$model <- model
    out$cotype <- cotype
    class(out) <- "summary.prodlim"
    out
    # }}}
}
