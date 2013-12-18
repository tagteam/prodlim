## Used to specify nearest neighborhoods for a continuous predictor.
## Similar as the function strata specifies strata for a discrete predictor. 


#' Dummy function
#' 
#' This is a special function used in the context of the conditional product
#' limit estimator. It is similar to the function strata and can be used in a
#' formula to force a stratified analysis in nearest neighborhoods defined by
#' the values of a continuous predictor.
#' 
#' 
#' @param x A continuous predictor variable
#' @return The argument is passed on as it is: NN(x)=x
#' @author Thomas Gerds
#' @seealso \code{\link{prodlim}}
#' @keywords nonparametric
#' @examples
#' 
#' dat <- SimSurv(20)
#' prodlim(Hist(time,status)~NN(X2),data=dat)
#'
#' @export
NN <- function(x)x
