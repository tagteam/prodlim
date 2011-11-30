resolveLinPred <- function(X,coef,transform,transName="f",verbose=TRUE){
  if (is.null(X) || is.null(coef)) {
    LP <- 0
  }
  else {
    NP <- NCOL(X)
    NC <- length(coef)
    stopifnot((length(coef)>0) && all(is.numeric(coef)))
    if (NP != length(coef)){
      if (length(coef)==1){
        if (verbose) warning("The regression coefficient ",coef," is used for all covariates.")
        coef <- rep(coef,NP)
      }
      else{
        stop("Number of covariates and number of regression coefficients differ.")
      }
    }
    LP <- colSums(coef * t(X))
  }
  LP
}
