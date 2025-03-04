### as.data.table.prodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (12:57) 
## Version: 
## Last-Updated: Mar  3 2025 (14:11) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Collect results of a fitted prodlim object in a data.table
##'
##' By default object contains results for all fitted time points and all strata.
##' Use arguments times and newdata of \code{\link{summary.prodlim}} to subset.
##' @title Turn prodlim objects into a \code{\link[data.table]{data.table}}
##' @name as.data.table.prodlim
##' @param x object obtained with function \code{\link{prodlim}}
##' @param keep.rownames Not used
##' @param ... passed to \code{\link{summary.prodlim}}
##' @return A data.table with the results of the prodlim object 
##' @seealso \code{\link{prodlim}}, \code{\link[data.table]{data.table}}
##' @examples
##' library(data.table)
##' set.seed(8)
##' d <- SimCompRisk(17)
##' fit <- prodlim(Hist(time,event)~X1,data=d)
##' as.data.table(fit)
##'
##'@export as.data.table.prodlim
##'@export 
##'@author Thomas A. Gerds <tag@@biostat.ku.dk>
as.data.table.prodlim <- function(x,keep.rownames = FALSE,...){
    requireNamespace("data.table")
    all_times = sort(unique(c(0,x$time)))
    all_X = x$X
    out = data.table::as.data.table(summary(x,times = all_times,newdata = all_X,...,format = "df"))
    data.table::setnames(out,sub("cuminc","absolute_risk",names(out)))
    out[]
}


######################################################################
### as.data.table.prodlim.R ends here
