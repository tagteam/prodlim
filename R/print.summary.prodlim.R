##' @S3method print summary.prodlim
##' @method print summary.prodlim
print.summary.prodlim <- function(x,...){
    model <- x$model
    cotype <- x$cotype
    sumtable <- x$table
    if (class(sumtable)=="matrix"){
        print(sumtable,quote=FALSE,...)
    }
    else{
        if (model=="survival"){
            if (cotype==1){
                print(sumtable,quote=FALSE,...)
            } else{
                print.listof(sumtable,quote=FALSE,...)
            }
        } else{
            if (model=="competing.risks"){
                for (cc in 1:length(sumtable)){
                    cat("\n\n----------> Cause: ",names(sumtable)[cc],"\n\n")
                    if (cotype==1){
                        print(sumtable[[cc]],quote=FALSE,...)
                    }
                    else{
                        print.listof(sumtable[[cc]],quote=FALSE,...)
                    }
                }
            }
        }
    }
}
