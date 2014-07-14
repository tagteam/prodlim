#' Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
#' 
#' Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
#' 
#' 
#' @param x Object of class \code{"prodlim"}.
#' @param q Quantiles. Vector of values between 0 and 1.
#' @param ... not used
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @keywords survival
#' @examples
#' 
#' set.seed(1)
#' d=SimSurv(30)
#' f=prodlim(Hist(time,status)~X1,data=d)
#' # median
#' quantile(f,.5)
#'
#' @method quantile prodlim
#' @S3method quantile prodlim
"quantile.prodlim" <- function(x,
                               q,
                               ...){
  ## require(stats)
  stopifnot(x$model=="survival")
  if (missing(q)) q <- c(1,.75,0.5,.25,0)
  q <- 1-q ## since this is a survival function
  sumx <- summary(x,newdata=x$X,times=x$time,showTime=TRUE,verbose=FALSE)
  getQ <- function(sum){
    out <- do.call("cbind",lapply(c("surv","lower","upper"),function(w){
      notna=is.na(sum[,w])
      xxx=sum[,w][!notna]
      ttt=sum[,"time"][!notna]
      found <- 2+sindex(jump.times=xxx,eval.times=q,comp="greater",strict=FALSE)
      inner <- c(as.vector(c(0,ttt)[found]))
      inner
    }))
    out <- data.frame(out)
    out <- cbind(q,out)
    names(out) <- c("q","quantile","lower","upper")
    out}
   if (sumx$cotype==1) out <- list("quantiles.survival"=getQ(sumx$table))
  else out <- lapply(sumx$table,getQ)
  class(out) <- "quantile.prodlim"
  out
}

