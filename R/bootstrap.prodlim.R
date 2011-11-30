bootstrap.prodlim <- function(object,times,B){
  stopifnot(!(object$covariate.type==1 & object$model=="survival"))
}
