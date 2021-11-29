#' Add point-wise confidence limits to the graphs of Kaplan-Meier and
#' Aalen-Johansen estimates.
#' 
#' This function is invoked and controlled by \code{plot.prodlim}.
#' 
#' This function should not be called directly. The arguments can be specified
#' as \code{Confint.arg} in the call to \code{plot.prodlim}.
#' 
#' @param ci A \code{data.table} with columns \code{time}, \code{lower} and \code{upper}. 
#' @param citype If \code{"shadow"} then confidence limits are drawn as colored
#' shadows.  Otherwise, dotted lines are used to show the upper and lower
#' confidence limits.
#' @param col the colour of the lines.
#' @param lty the line type of the lines.
#' @param lwd the line thickness of the lines.
#' @param density For \code{citype="shadow"}, the density of the shade. Default
#' is 55 percent.
#' @param \dots Further arguments that are passed to the function
#' \code{segments} if \code{type=="bars"} and to \code{lines} else.
#' @return Nil
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot.prodlim}}, \code{\link{atRisk}},
#' \code{\link{markTime}}
#' @keywords survival
#' @export
confInt <- function(ci,
                    citype,
                    col,
                    lty,
                    lwd,
                    density=55,
                    ...){
    lower <- upper <- time <- NULL
    nlines <- length(ci)
    nix <- lapply(1:nlines,function(i){
        # ----------remove confidence limits before the first event----------
        CI <- ci[[i]][upper-lower<1]
        if (NROW(ci[[i]])>0){
            switch(citype,
                   "bars"={
                       nix <- CI[,{segments(x0=time,x1=time,y0=lower,y1=upper,lwd=lwd[i],col=col[i],lty=lty[i],...);NULL}]
                   },
                   "shadow"={
                       nix <- CI[,{
                           cc <- dimColor(col[i],density=density)
                           ttt <- time
                           nt <- length(ttt)
                           ttt <- c(ttt,ttt)
                           uuu <- c(0,upper[-nt],upper)
                           lll <- c(0,lower[-nt],lower)
                           neworder <- order(ttt)
                           uuu <- uuu[neworder]
                           lll <- lll[neworder]
                           ttt <- sort(ttt)
                           polygon(x=c(ttt,rev(ttt)),
                                   y=c(lll,rev(uuu)),col=cc,border=NA)
                           NULL}]
                   },{
                       nix <- CI[,{
                           lines(x=time,lower,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
                           lines(x=time,upper,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
                           NULL
                       }]
                   })
        }
    })
}
