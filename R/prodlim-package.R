### prodlim-package.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Apr 24 2015 (09:08) 
## Version: 
## last-updated: Jun 24 2024 (11:58) 
##           By: Thomas Alexander Gerds
##     Update #: 30
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Kaplan-Meier and Aalen-Johansen method
#'
#' Fast and user friendly implementation
#' of nonparametric estimators for censored
#' event history (survival) analysis.
#' @keywords internal
# "_PACKAGE"
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
