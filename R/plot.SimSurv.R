plot.SimSurv <- function(x,...){
  form <- attr(x,"formula")
  fit <- prodlim(form,data=x)
  plot(fit,...)
}
