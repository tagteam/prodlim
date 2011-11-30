find.baseline <- function(x=.5,
                          setting,
                          verbose=FALSE){
  N <- setting$N
  f <- function(y){
    setting$cens.baseline <- y
    ncens <- sum(do.call("SimSurv",replace(setting,"verbose",verbose))$status==0)
    x-ncens/N
  }
  base.cens <- uniroot(f,c(exp(-50),1000000),tol=.0000001,maxiter=100)$root
  new.setting <- setting
  new.setting$cens.baseline <- base.cens
  do.call("SimSurv",replace(new.setting,"verbose",TRUE))
  new.setting
}
