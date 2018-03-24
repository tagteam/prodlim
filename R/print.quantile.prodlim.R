##' @export 
print.quantile.prodlim <- function(x,digits=2,na.val="--",...){
    if (attr(x,"reverse")==FALSE)
        cat("Quantiles of the event time distribution based on the ",
            ifelse(x$type=="surv","Kaplan-Meier","Aalen-Johansen"),
            "method.")
    else
        cat("Quantiles of the potential follow up time distribution based on the Kaplan-Meier method",
            "\napplied to the censored times reversing the roles of event status and censored.")
    thisfmt <- paste0("%1.",digits[[1]],"f")
    printx <- function(u){
        ifelse(is.na(u),na.val,sprintf(thisfmt,u))
    }
    cat("\n")
    lapply(1:length(x),function(i){
        tab <- x[[i]]
        cat("\nTable of quantiles and corresponding confidence limits:\n")
        if (attr(x,"cotype")!=1) {
            if(length(names(x)[i])){
                cat("\n",names(x)[i],"\n\n")
            }
        }
        print(tab,digits=digits)
        ## if(0.5 %in% tab$q){
        ## cat("Median time (lower.95; upper.95): ",printx(tab[tab$q==0.5,"quantile"])," (",printx(tab[tab$q==0.5,"lower"]),";",printx(tab[tab$q==0.5,"upper"]),")","\n",sep="")
        ## }
        if(all(c(0.25,0.5,0.75) %in% tab$q)){
            if (attr(x,"model")=="survival")
                cat("Median time (IQR):",printx(tab[tab$q==0.5,"quantile"])," (",printx(tab[tab$q==0.75,"quantile"]),";",printx(tab[tab$q==0.25,"quantile"]),")","\n",sep="")
            else
                cat("Median time (IQR):",printx(tab[tab$q==0.5,"quantile"])," (",printx(tab[tab$q==0.25,"quantile"]),";",printx(tab[tab$q==0.75,"quantile"]),")","\n",sep="")
        }
    })
    invisible(x)
}
