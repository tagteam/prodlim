### prodlim-package.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Apr 24 2015 (09:08) 
## Version: 
## last-updated: Mar  6 2017 (12:32) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Functions for estimating probabilities from right censored data
#'
#' @docType package
#' @name prodlim
#' @useDynLib prodlim, .registration=TRUE
#' @importFrom survival survdiff Surv cluster
#' @importFrom stats quantile 
#' @import lava
#' @importFrom Rcpp sourceCpp
## --> importFrom KernSmooth dpik
#' @importFrom graphics abline axis lines mtext par plot points polygon rect segments strheight strwidth text
#' @importFrom stats .getXlevels delete.response drop.terms formula get_all_vars median model.frame model.matrix model.response na.omit pchisq predict qnorm reformulate terms update update.formula
NULL


#----------------------------------------------------------------------
### prodlim-package.R ends here
