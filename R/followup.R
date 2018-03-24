### followup.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 22 2015 (10:29) 
## Version: 
## last-updated: Sep 23 2017 (14:00) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
followup <- function(formula,data,...){
    G <- prodlim(formula,data,reverse=TRUE)
    quantile(G,...)
}


#----------------------------------------------------------------------
### followup.R ends here
