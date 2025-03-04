### as.data.frame.prodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (12:57) 
## Version: 
## Last-Updated: Mar  3 2025 (14:13) 
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
##' 
##'@export as.data.frame.prodlim
##'@export 
##'@author Thomas A. Gerds <tag@@biostat.ku.dk>
as.data.frame.prodlim <- function(x,...){
    all_times = sort(unique(c(0,x$time)))
    all_X = x$X
    out = summary(x,times = all_times,newdata = all_X,...,format = "df")
    class(out) = "data.frame"
    names(out) = sub("cuminc","absolute_risk",names(out))
    out[]
}


######################################################################
### as.data.frame.prodlim.R ends here
