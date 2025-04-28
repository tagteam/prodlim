### as.data.frame.prodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (12:57) 
## Version: 
## Last-Updated: Mar  5 2025 (16:14) 
##           By: Thomas Alexander Gerds
##     Update #: 41
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Collect results of a fitted prodlim object in a data.frame
##'
##' By default object contains results for all fitted time points and all strata.
##' Use arguments times and newdata of \code{\link{summary.prodlim}} to subset.
##' @title Turn prodlim objects into a data.frame
##' @name as.data.frame.prodlim
##' @param x object obtained with function \code{\link{prodlim}}
##' @param ... passed to \code{\link{summary.prodlim}}
##' @return A data.table with the results of the prodlim object 
##' @seealso \code{\link{prodlim}}
##' @examples
##' set.seed(8)
##' d <- SimCompRisk(17)
##' fit <- prodlim(Hist(time,event)~X1,data=d)
##' as.data.frame.prodlim(fit)
##' as.data.frame.prodlim(fit)
##' 
##'@export as.data.frame.prodlim
##'@export 
##'@author Thomas A. Gerds <tag@@biostat.ku.dk>
as.data.frame.prodlim <- function(x,...){
    args <- list(...)
    if ("newdata"%in%names(args)) X <- args$newdata else X <- x$X 
    if ("times"%in%names(args)) ttt <- args$times else ttt <- sort(unique(c(0,x$time)))
    if ("cause"%in%names(args)) ccc <- args$cause else ccc <- attr(x$model.response,"states")
    args <- c(list(object = x,
                   times = ttt,
                   newdata = X,
                   cause = ccc,
                   format = "df"),args)
    out <- do.call("summary.prodlim",args[!duplicated(names(args))])
    class(out) = "data.frame"
    names(out) = sub("cuminc","absolute_risk",names(out))
    out[]
}


######################################################################
### as.data.frame.prodlim.R ends here
