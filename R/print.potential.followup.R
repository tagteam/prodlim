### print.potential.followup.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (08:23) 
## Version: 
## Last-Updated: Mar  3 2025 (14:12) 
##           By: Thomas Alexander Gerds
##     Update #: 35
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' @export 
print.potential.followup <- function(x,digits=2,na.val="--",...){
    Q1 = Q3 = NULL
    cat("Median and inter quartile range (IQR) of the potential follow up time based on the reverse Kaplan-Meier method.\n")
    fmt = paste0("%1.",digits,"f")
    class(x) = c("data.table","data.frame")
    sx = copy(x)
    set(sx,j = "Median (IQR)",value = gsub("NA","--",sprintf(fmt = paste0(fmt," (",fmt,";",fmt,")"),x$median,x$Q3,x$Q1)))
    psx = copy(sx)
    psx[,median := NULL]
    psx[,Q1 := NULL]
    psx[,Q3 := NULL]
    print(psx)
    invisible(sx)
}


######################################################################
### print.potential.followup.R ends here
