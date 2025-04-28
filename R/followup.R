### followup.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 22 2015 (10:29) 
## Version: 
## last-updated: Apr 25 2025 (08:17) 
##           By: Thomas Alexander Gerds
##     Update #: 37
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' The reverse Kaplan-Meier method estimates the median potential followup time. 
##'
##' This is merely a wrapper for \code{prodlim} with argument reverse. 
##' @title Estimation of the median potential followup time.
##' @param formula A formula whose left hand side is a \code{Hist} or a 
##' \code{\link[survival]{Surv}} object specifying the event time
##' and the event type where by default 0=censored, 1=event, 2=competing risk (if any).
##' Use \code{cens.code} to change the value for censored.
##' @param cens.code Value of the event 
##' @param data A data.frame in which all the variables of
##' \code{formula} can be interpreted.
##' @param ... Arguments passed to \code{prodlim.quantile}. 
##' @return
##' The estimated median potential followup time with inter quartile ranges. 
##' @seealso \code{\link{prodlim}}
##' @examples
##' set.seed(8)
##' d <- SimCompRisk(117)
##'
##' # overall
##' followup(Hist(time,event)~1,data=d)
##' 
##' # in strata defined by variable X1 
##' followup(Hist(time,event)~X1,data=d)
##' @references Michael Schemper and Terry L. Smith. A note on quantifying follow-up in studies of failure time. Controlled Clinical Trials, 17(4):343--346, 1996.
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
followup <- function(formula,cens.code = 0,data,...){
    G <- prodlim(formula,data,reverse=TRUE)
    x <- quantile(G,q = c(0.25,0.5,0.75),...)
    x <- data.table::copy(x)
    class(x) = c("data.table","data.frame")
    data.table::set(x,j = "lower",value = NULL)
    data.table::set(x,j = "upper",value = NULL)
    class(x) <- c("potential.followup","data.table","data.frame")
    byvars <- c(G$discrete.predictors,G$continuous.predictors)
    if (length(byvars) == 0) {
        byvars <- "X"
        x <- cbind(X = "Overall",x)
    }
    ff <- formula(paste0(paste0(byvars,collapse = "+"),"~q"))
    x <- data.table::dcast(data = x,formula = ff,value.var = "quantile")
    setnames(x,c("0.25","0.5","0.75"),c("Q1","median","Q3"))
    attr(x,"covariates") <- byvars
    class(x) = c("potential.followup","data.table","data.frame")
    x
}
#----------------------------------------------------------------------
### followup.R ends here
