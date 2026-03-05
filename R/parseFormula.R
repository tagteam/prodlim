### parseFormula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: mar  4 2026 (11:28) 
## Version: 
## Last-Updated: mar  4 2026 (14:33) 
##           By: Thomas Alexander Gerds
##     Update #: 23
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Parse a formula to extract terms and handle special functions
##'
##' This function parses a model formula, extracts terms, and processes specified 'special' functions and their arguments.
##' It allows defining default values for special functions and specifying aliases for those names.
##'
##' @title Parse Formula
##' @param formula A model formula (e.g., Y ~ X + Z or Y ~ fun(X, param = value)).
##' @param data An optional data frame containing the variables referenced in the formula.
##' @param specials A character vector of special function names that should be specifically handled.
##' @param special.default.values A named list of named lists specifying default argument values for each special.
##' @param alias.names A named list mapping special names to their alias(es) (e.g., `list("special" = "alias")`).
##' @return A list where each element corresponds to a special function in the formula and contains the parsed terms and their arguments.
##' @seealso [stats::terms()], [strip.terms()]
##' @examples
##' # Basic usage with specials
##' parseFormula(Y ~ fun(X, a = 1), specials = "fun")
##'
##' # With defaults and aliases
##' parseFormula(Y ~ gfun(X, a = 7), 
##'              specials = "fun", 
##'              alias.names = list("fun" = "gfun"),
##'              special.default.values = list("fun" = list("a" = 0)))
##'
##' # Complex use case with multiple specials
##' parseFormula(Surv(time, status) ~ age + strata(sex, test = 1) + prop(albumin),
##'              specials = c("prop", "strata"), 
##'              special.default.values = list("prop" = list(power = 0), "strata" = list(test = 1)))
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @return A named list with an element for each special that contains a named list where each element
##'         corresponds to a formula term that is affected by the special that contains
##'         the value(s) of the additional argument(s) for this term.
##' @seealso model.design strip.terms
##' @examples
##' parseFormula(Y~fun(X,a=1),specials="fun",special.default.values=list("fun"=list("a"=0)))
##' parseFormula(Y~fun(X),specials="fun",special.default.values=list("fun"=list("a"=0)))
##' parseFormula(Y~gfun(X,a=7),specials="fun",
##'              alias.names=list("fun"="gfun"),
##'              special.default.values=list("fun"=list("a"=0)))
##' parseFormula(formula = Surv(time,status)~age
##'                             +const(factor(edema))
##'                             +strata(sex,test=0)
##'                             +prop(alb)+prop(bili,power=1)
##'                             +tp(albumin),
##'                   specials = c("prop","timevar","strata","tp","const"),
##'                   special.default.values =list("prop" = list(power = 0),
##'                                                "strata" = list(test = 1)))
##' parseFormula(formula = Surv(time,status)~age
##'                             +const(factor(edema))
##'                             +strata(sex,test=0)
##'                             +prop(alb)+prop(bili,power=1)
##'                             +tp(albumin),
##'                   specials = c("prop","timevar","strata","tp"),
##'                   special.default.values =list("prop" = list(power = 0),
##'                                                "strata" = list(test = 1)),
##'                   alias.names=list("prop"="const"))
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
parseFormula <- function(formula,
                         data = NULL,
                         specials,
                         special.default.values = NULL,
                         alias.names = NULL){
    Terms <- stats::terms(x=formula,specials=specials,data=data)
    Terms <- strip.terms(Terms,
                         specials=specials,
                         arguments=special.default.values,
                         alias.names=alias.names,
                         unspecials = NULL)
    result <- list()
    for (special in specials) {
        strippedArgs <- attr(Terms, "stripped.arguments")[[special]]
        if (is.null(strippedArgs)) {
            result[[special]] <- NULL
        } else {
            parsedElements <- list()
            for (element in names(strippedArgs)) {
                elementArgs <- strippedArgs[[element]]
                if (is.null(elementArgs)) {
                    parsedElements[[element]] <- NULL
                } else {
                    parsedElements[[element]] <- lapply(elementArgs, function(arg) {
                        if (!is.na(as.numeric(arg))) as.numeric(arg) else arg
                    })
                }
            }
            result[[special]] <- parsedElements
        }
    }
    result
}



######################################################################
### parseFormula.R ends here
