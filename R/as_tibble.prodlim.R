### as_tibble.prodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (12:57) 
## Version: 
## Last-Updated: Mar  3 2025 (14:12) 
##           By: Thomas Alexander Gerds
##     Update #: 30
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
##' as_tibble(fit)
##'
##'@export as_tibble.prodlim
##'@export 
##'@author Thomas A. Gerds <tag@@biostat.ku.dk>
as_tibble.prodlim <- function(x,...){
    requireNamespace("tibble")
    all_times = sort(unique(c(0,x$time)))
    all_X = x$X
    out = summary(x,times = all_times,newdata = all_X,...,format = "df")
    class(out) = "data.frame"
    names(out) = sub("cuminc","absolute_risk",names(out))
    out = tibble::as_tibble(out)
    out[]
}


######################################################################
### as_tibble.prodlim.R ends here

