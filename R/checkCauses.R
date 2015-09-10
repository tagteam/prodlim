### checkCauses.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 10 2015 (11:56) 
## Version: 
## last-updated: Sep 10 2015 (11:59) 
##           By: Thomas Alexander Gerds
##     Update #: 2
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
    fitted.causes <- attributes(object$model.response)$states
    ## stopifnot(length(fitted.causes)==length(object$n.event))
    if (!is.numeric(cause)){
        Found <- match(as.character(cause),fitted.causes,nomatch=0)
        if (any(Found==0)) stop("Cannot find competing cause(s) ", as.character(cause)[Found==0], "in fitted object.")
        return(cause)
    }else{
         if (length(fitted.causes)<max(cause))
             stop(paste0("Object has fitted ",length(fitted.causes)," competing causes. So, there is no ",max(cause),"th cause."))
         return(fitted.causes[cause])
     }
}


#----------------------------------------------------------------------
### checkCauses.R ends here
