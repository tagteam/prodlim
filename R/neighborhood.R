"neighborhood" <- function(x,bandwidth=NULL,kernel="box"){
  N <- length(x)
  if (N<2) stop("Not enough observations for kernel smoothing.")
  
  orderx <- order(x)
  values <- sort(unique(x))
  NU <- length(values)
  workx <- factor(x,labels=1:NU)
  tabu <- tabulate(workx)
  cumtabu <- cumsum(tabu)
  cumtabx <- rep(cumtabu,tabu)
  tabx <- rep(tabu,tabu)
  if (!length(bandwidth)){ ## need a bandwidth (dpik is from KernSmooth)
    ## require(KernSmooth)
    bandwidth <- dpik(cumtabx/N,kernel="box")
  }
  else
    if (bandwidth=="smooth") bandwidth <- N^{-1/4}
  radius <- floor(bandwidth*N)
  
  nbh <- .C("neighborhood",first=integer(NU),size=integer(NU),as.integer(cumtabu),as.integer(cumtabx),as.integer(tabx),as.integer(radius),as.integer(NU),as.integer(N),PACKAGE="prodlim")
  nall <- sum(nbh$size)
  neighbors <- .C("neighbors",first=nbh$first,size=nbh$size,as.integer(orderx),neighbors=integer(nall),as.integer(NU),PACKAGE="prodlim")$neighbors
  
  out <- list(values=values,
              first.nbh=nbh$first,
              size.nbh=nbh$size,
              neighbors=neighbors,
              bandwidth=bandwidth,
              kernel=kernel,
              nu=NU,
              n=N)
    
  class(out) <- "neighborhood"
  out
}
