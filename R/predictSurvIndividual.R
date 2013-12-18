predictSurvIndividual <- function(object,
                                  lag=1){
  obs.times <- as.numeric(object$model.response[,1])
  if (object$covariate.type==1){
    locOBS <- match(obs.times,object$time,nomatch=FALSE)
    if (any(locOBS==FALSE)) stop("Can't locate all individual observation times" )
    psurv <- c(rep(1,lag),object$surv)[locOBS]}
  else{
    N <- length(obs.times)
    findex <- row.match(object$model.matrix,object$X)
    psurv <- .C("predict_individual_survival",pred=double(N),as.double(object$surv),as.double(object$time),as.double(obs.times),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(N),as.integer(lag),NAOK=FALSE,PACKAGE="prodlim")$pred}
  psurv
}
