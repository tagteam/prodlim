##' @export 
print.quantile.prodlim <- function(x,byvars,digits=2,na.val="--",...){
    if (attr(x,"reverse")==FALSE)
        cat("Quantiles of the event time distribution based on the ",
            ifelse(attr(x,"model")=="survival","Kaplan-Meier","Aalen-Johansen"),
            " method.",sep="")
    else
        cat("Quantiles of the potential follow up time distribution based on the Kaplan-Meier method",
            "\napplied to the censored times reversing the roles of event status and censored.")
    cat("\n")
    cat("\nTable of quantiles and corresponding confidence limits:\n")
    tab <- data.table::data.table(do.call("data.frame",x))
    print(tab,digits=digits)
    if(all(c(0.25,0.5,0.75) %in% tab$q)){
        set(tab,j = "quantile",value = sprintf("%1.2f",tab[["quantile"]]))
        byvars = attr(x,"covariates")
        cat("\n")
        if (attr(x,"model")=="survival"){
            tab[,{
                if (length(byvars)>0) cat(paste(byvars,"=",.SD[1,byvars,with = FALSE],collapse = ", "),":")
                cat("Median time (IQR):",quantile[q==0.5]," (",quantile[q==0.75],";",quantile[q==0.25],")","\n",sep="")},by = byvars,.SDcols = c(byvars,"q","quantile")
                ]
        }
        else{
            tab[,{
                if (length(byvars)>0) cat(paste(byvars,"=",.SD[1,byvars,with = FALSE],collapse = ", "),":")
                cat("Median time (IQR):",quantile[q==0.5]," (",quantile[q==0.25],";",quantile[q==0.75],")","\n",sep="")},by = byvars,.SDcols = c(byvars,"q","quantile")
                ]
        }
    }
    invisible(x)
}
