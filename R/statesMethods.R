states <- function(object,...){
  UseMethod("states",object)
}

states.Hist <- function(object){
  attr(object,"states")
}


states.prodlim <- function(object){
  attr(object$model.response,"states")
}

