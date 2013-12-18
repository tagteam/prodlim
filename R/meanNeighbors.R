meanNeighbors <- function(x,y,...){
  nnn=neighbors(x,y,...)
  out <- data.frame(x=nnn$nbh$values,
                    y=sapply(nnn$list,mean))
  names(out) <- c("uniqueX","averageY")
  out
}
