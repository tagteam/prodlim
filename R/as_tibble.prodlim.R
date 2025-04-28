### as_tibble.prodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (12:57) 
## Version: 
## Last-Updated: Mar  5 2025 (15:42) 
##           By: Thomas Alexander Gerds
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Collect results of a fitted prodlim object in a tibble
##'
##' By default object contains results for all fitted time points and all strata.
##' Use arguments times and newdata of \code{\link{summary.prodlim}} to subset.
##' @title Turn prodlim objects into a \code{\link[tibble]{tibble}}
##' @name as_tibble.prodlim
##' @param x object obtained with function \code{\link{prodlim}}
##' @param ... passed to \code{\link{summary.prodlim}}
##' @return A data.table with the results of the prodlim object 
##' @seealso \code{\link{prodlim}}, \code{\link[tibble]{tibble}}
##' @examples
##' library(tibble)
##' set.seed(8)
##' d <- SimCompRisk(17)
##' fit <- prodlim(Hist(time,event)~X1,data=d)
##' tibble::as_tibble(fit)
##'
##'@export as_tibble.prodlim
##'@export 
##'@author Thomas A. Gerds <tag@@biostat.ku.dk>
as_tibble.prodlim <- function(x,...){
    requireNamespace("tibble")
    tibble::as_tibble(as.data.frame.prodlim(x,...))
}


######################################################################
### as_tibble.prodlim.R ends here

