##' Extract the states of a multi-state model
##'
##' Applying this function to the fit of prodlim means to apply
##' it to \code{fit$model.response}.
##' @title States of a multi-state model
##' @param object
##' @param ...
##' @return A character vector with the states of the model.
##' @author Thomas A. Gerds
#' @export
getStates <- function(object,...){
  UseMethod("getStates",object)
}
#' @method getStates Hist
#' @S3method getStates Hist
getStates.Hist <- function(object,...){
  attr(object,"states")
}

#' @method getStates prodlim
#' @S3method getStates prodlim
getStates.prodlim <- function(object){
  attr(object$model.response,"states")
}

