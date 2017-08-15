### checkCauses.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 10 2015 (11:56) 
## Version: 
## last-updated: Aug 15 2017 (06:43) 
##           By: Thomas Alexander Gerds
##     Update #: 10
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
    if (!is.character(cause)) cause <- as.character(cause)
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

#----------------------------------------------------------------------
### checkCauses.R ends here
