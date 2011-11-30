shiftPosition <- function(A,dir,offset){
  ## if (dir %*% c(0,1) ==0) # vertical
  A + sign(dir) * c(offset,offset)
}

## shiftPosition2 <- function(A,B,offset){
##   out <- list(x1=NULL,y1=NULL,x2=NULL,y2=NULL)
##   if (A[1]==B[1]){# A vertically above or below B
##     out$x1 <- out$x2 <- A[1]
##     if (A[2]<B[2]) { # A below B
##       out$y1 <- A[2] + offset
##       out$y2 <- B[2] - offset
##     }
##     else{ # A above B
##       out$y1 <- A[2] - offset
##       out$y2 <- B[2] + offset
##     }
##   }
##   else{
##     thisform <- function(x,y){
##       (y[4]-y[2])*x/(y[3]-y[1]) - (y[4]-y[2])*y[1]/(y[3] -y[1]) + y[2]
##     }
##     out$x1 <- A[1] + offset
##     out$x2 <- B[1] - offset
##     out$y1 <- thisform(x=out$x1,y=c(A,B))
##     out$y2 <- thisform(x=out$x2,y=c(A,B))
##   }
##   out
## }
