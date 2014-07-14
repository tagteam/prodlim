##' Extract event history data and design matrix including specials from call
##'
##' Obtain a list with the data used for event history regression analysis. This
##' function cannot be used directly on the user level but inside a function
##' to prepare data for survival analysis. 
##' @title Event history frame
##' @param call
##' @param specials 
##' @return 
##' @seealso model.frame model.design 
##' @examples
##' ## event time no competing risks
##' set.seed(13)
##' dsurv <- SimSurv(4)
##' dsurv$X1 <- factor(dsurv$X1)
##' formsurv <- Hist(time,status)~afun(X1)+X2
##' prop <- afun <- function(x)x
##' areg <- function(formula,data=parent.frame()){
##'    thecall <- match.call()
##'    EventHistory.frame(formula=formula,data=data,specials=c("prop","afun"))
##' }
##' areg(formsurv,dsurv)
##' time = 1:2
##' status = c(0,1)
##' X2=X1=c(4,5)
##' areg(formsurv)
##'
##' formsurv2 <- Hist(time,status)~prop(X1)
##' areg(formsurv2)
##' ff <- list(formsurv,formsurv2)
##' lapply(ff,function(f){areg(f,dsurv)})
##' 
##' set.seed(170)
##' dcr <- SimCompRisk(10)
##' dcr$entry <- dcr$time-dcr$time/2
##' formcr <- Hist(time,event)~X1+prop(X2)
##' formcr <- Hist(time,entry=entry,event=event)~X1+prop(X2)
##' areg(formcr,dcr)
##'
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
EventHistory.frame <- function(formula,data,specials,specials.to.factor="strata"){
    ## define specials
    ## for (s in specials)
        ## assign(s,function(x)x)
    Terms <- terms(x=formula, specials=specials, data = data)
    if (sum(sapply(c("Surv","Hist"),function(x) length(grep(x,Terms[[2]])>0)))==0)
        stop("Expected a 'Surv'-object or a 'Hist'-object on the left hand side of the formula.")
    m <- model.frame(formula=Terms,data=data)
    ## m$formula <- Terms
    ## m[[1]] <- as.name("model.frame")
    ## m <- eval(m, parent.frame())
    design <- model.design(data=m,
                           max.order=1,
                           drop.intercept=TRUE,
                           specials.to.factor=specials.to.factor)
    X <- design$X
    strata <- design$strata[[1]]
    id <- design$cluster
    Y <- model.extract(m, "response")
    # Fix for those who use `Surv' instead of `Hist' 
    if (match("Surv",class(Y),nomatch=0)!=0){
        attr(Y,"model") <- "survival"
        attr(Y,"cens.type") <- "rightCensored"
        attr(Y,"entry.type") <- ifelse(ncol(Y)==2,"","leftTruncated")
        if (attr(Y,"entry.type")=="leftTruncated")
            colnames(Y) <- c("entry","time","status")
    }
    delayed <- attr(Y,"entry.type")=="leftTruncated"
    class(Y) <- "matrix"
    eth <- lapply(data.frame(Y),function(x)x)
    if (attr(Y,"model")=="competing.risk")
        eth$event <- getEvent(Y,mode="factor")
    out <- list(eth,design)
    out
}
