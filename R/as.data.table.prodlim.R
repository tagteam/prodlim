### as.data.table.prodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (12:57) 
## Version: 
## Last-Updated: Mar  5 2025 (11:58) 
##           By: Thomas Alexander Gerds
##     Update #: 23
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
    data.table::as.data.table(as.data.frame.prodlim(x,...))
}


######################################################################
### as.data.table.prodlim.R ends here
