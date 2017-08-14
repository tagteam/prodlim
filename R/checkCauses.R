### checkCauses.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 10 2015 (11:56) 
## Version: 
## last-updated: Jul  7 2017 (11:05) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
checkCauses <- function(cause,object){
    cause <- unique(cause)
    fitted.causes <- attr(object$model.response,which="states")
    if (!(all(cause %in% fitted.causes))){
        stop(paste0("Cannot find all requested cause(s) in attr(object$model.response,'states')\n\n",
                     "Requested cause(s): ",
                     paste0(cause,collapse=", "),
                     "\n Available causes: ",
                     paste(fitted.causes,collapse=", "),"\n"))
    }
    cause
}
## stopifnot(length(fitted.causes)==length(object$n.event))
## if (!is.numeric(cause)){
## Found <- match(as.character(cause),fitted.causes,nomatch=0)
## if (any(Found==0)) stop("Cannot find competing cause(s) ", as.character(cause)[Found==0], "in fitted object.")
## return(cause)
## }else{
## if (length(fitted.causes)<max(cause))
## stop(paste0("Object has fitted ",length(fitted.causes)," competing causes. So, there is no cause number: ",max(cause)))
## return(fitted.causes[cause])
## }
## }


#----------------------------------------------------------------------
### checkCauses.R ends here
