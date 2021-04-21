##' @export 
print.summary.prodlim <- function(x,digits,...){
    if (!is.null(attr(x,"model"))){
        if (missing(digits))
            digits <- ifelse(attr(x,"percent"),1,3)
        print.data.frame(x,digits=digits,...)
    } else{ # old behaviour from before May 2021
        model <- x$model
        cotype <- x$cotype
        sumtable <- x$table
        if (missing(digits))
            digits <- ifelse(x$percent,1,3)
        if (model=="survival"){
            if (cotype==1){
                print(sumtable,digits=digits,quote=FALSE,...)
            } else{
                print.listof(sumtable,digits=digits,quote=FALSE,...)
            }
        } else{
            if (model=="competing.risks"){
                for (cc in 1:length(sumtable)){
                    cat("\n\n----------> Cause: ",names(sumtable)[cc],"\n\n")
                    if (cotype==1){
                        print(sumtable[[cc]],digits=digits,quote=FALSE,...)
                    }
                    else{
                        print.listof(sumtable[[cc]],digits=digits,quote=FALSE,...)
                    }
                }
            }
        }
    }
}
