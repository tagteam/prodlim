### prodlim-package.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Apr 24 2015 (09:08) 
## Version: 
## last-updated: mar  6 2026 (08:49) 
##           By: Thomas Alexander Gerds
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
utils::globalVariables(c(
  "xmin","xmax","ymin","ymax","cause",".gg_kind",".gg_groupid",".gg_colour",".gg_fill"
))
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
