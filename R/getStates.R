##' Extract states from prodlim object
##'
##' Extract the two states of a simple survival model (initial and event) or the multiple states from a competing risks or a multi-state model.
##' @title Extract states 
##' @param object Object of class \code{"Hist"}.
##' @return states 
##' @author Thomas Alexander Gerds
##' @export
getStates <- function(object){
  attr(object,"states")
}
