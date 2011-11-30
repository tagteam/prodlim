PercentAxis <- function(x,at,...){
  axis(x,at=at,labels=paste(100*at,"%"),...)
}
