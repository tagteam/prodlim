#' Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
#' 
#' Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
#' 
#' 
#' @param x Object of class \code{"prodlim"}.
#' @param q Quantiles. Vector of values between 0 and 1.
#' @param cause For competing risks the cause of interest.
#' @param ... not used
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @keywords survival
#' @examples
#' library(lava)
#' set.seed(1)
#' d=SimSurv(30)
#'
#' # Quantiles of the potential followup time
#' g=prodlim(Hist(time,status)~1,data=d,reverse=TRUE)
#' quantile(g)
#' 
#' # survival time
#' f=prodlim(Hist(time,status)~1,data=d)
#' f1=prodlim(Hist(time,status)~X1,data=d)
#' # default: median and IQR
#' quantile(f)
#' quantile(f1)
#' # median alone
#' quantile(f,.5)
#' quantile(f1,.5)
#'
#' # competing risks
#' set.seed(3)
#' dd = SimCompRisk(30)
#' ff=prodlim(Hist(time,event)~1,data=dd)
#' ff1=prodlim(Hist(time,event)~X1,data=dd)
#' ## default: median and IQR
#' quantile(ff)
#' quantile(ff1)
#' 
#' print(quantile(ff1),na.val="NA")
#' print(quantile(ff1),na.val="Not reached")
#' 
#' @export quantile.prodlim
#' @export 
"quantile.prodlim" <- function(x,
                               q,
                               cause=1,
                               ...){
    time <- surv <- lower <- upper <- cuminc <- NULL
    get.quantiles <- function(time,x,lower,upper,model="survival"){
        out <- do.call("cbind",lapply(list(x,lower,upper),function(sumw){
            notna= is.na(sumw) | sumw==0 | sumw ==1
            if (all(notna)) return(NA)
            xxx=as.numeric(sumw[!notna])
            ttt=as.numeric(time[!notna])
            found <- 2+sindex(jump.times=xxx,eval.times=q,comp=ifelse(model=="survival","greater","smaller"),strict=FALSE)
            inner <- c(as.vector(c(0,ttt)[found]))
            inner
        }))
        out <- data.frame(out)
        out <- cbind(q,out)
        if (model=="survival") {
            names(out) <- c("q","quantile","lower","upper")
        }else{
            names(out) <- c("q","quantile","upper","lower")
            out <- out[,c("q","quantile","lower","upper")]
        }
        out
    }
    etype <- attr(x$model.response,"entry.type")
    if (!is.null(etype) && etype=="leftTruncated")
        stop("Don't know how to compute quantiles with delayed entry (left-truncation).")
    if(x$model=="survival"){
        if (missing(q)) q <- c(1,.75,0.5,.25,0)
        q <- 1-q ## since this is a survival function
        sumx <- summary(x,newdata=x$X,times=x$time,verbose=FALSE)
        if (attr(sumx,"cotype")==1) {
            out <- sumx[,get.quantiles(time=time,x=surv,lower=lower,upper=upper)]
        } else{
            out <- sumx[,get.quantiles(time=time,x=surv,lower=lower,upper=upper,model="survival"),by=key(sumx)]
        }
    } else{
        ## absolute risks, cumulative incidence, competing risks
        if (missing(q)) q <- c(0,0.25,0.5,0.75,1)
        sumx <- summary(x,newdata=x$X,times=x$time,verbose=FALSE,cause=cause)
        out <- sumx[,get.quantiles(time=time,x=cuminc,lower=lower,upper=upper,model = "cuminc"),by = key(sumx)]
        out
    }
    attr(out,"model") <- x$model
    attr(out,"covariates") <- key(sumx)
    attr(out,"reverse") <- x$reverse
    attr(out,"cotype") <- attr(sumx,"cotype")
    class(out) <- "quantile.prodlim"
    out
}
      

