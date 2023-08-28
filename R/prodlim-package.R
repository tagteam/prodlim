### prodlim-package.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Apr 24 2015 (09:08) 
## Version: 
## last-updated: Aug 28 2023 (09:18) 
##           By: Thomas Alexander Gerds
##     Update #: 26
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Fast and user friendly implementation
#' of nonparametric estimators for censored
#' event history (survival) analysis.
#' Kaplan-Meier and Aalen-Johansen method.
#'
#' @keywords internal
# "_PACKAGE"
#' @docType package
#' @name prodlim-package
#' @aliases prodlim-package
#' @useDynLib prodlim, .registration=TRUE
#' @importFrom survival survdiff Surv cluster
#' @import data.table
#' @importFrom stats quantile
#' @importFrom grDevices rainbow 
#' @import lava
#' @importFrom Rcpp sourceCpp
## --> importFrom KernSmooth dpik
#' @importFrom graphics grconvertX grconvertY abline axis lines mtext par plot points polygon rect segments strheight strwidth text
#' @importFrom stats .getXlevels delete.response drop.terms formula get_all_vars median model.frame model.matrix model.response na.omit pchisq predict qnorm reformulate terms update update.formula
NULL


#----------------------------------------------------------------------
### prodlim-package.R ends here
