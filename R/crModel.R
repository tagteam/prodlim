#' Competing risks model for simulation
#' 
#' Create a competing risks model with to causes to simulate a right censored event time data without
#' covariates
#' 
#' This function requires the \code{lava} package.
#' @title Competing risks model for simulation
#' @return A structural equation model initialized with four variables: the
#' latent event times of two causes, the latent right censored time, and the observed
#' right censored event time.
#' @author Thomas A. Gerds
#' @export
crModel <- function(){
    # require(lava)
    crm <- lvm()
    distribution(crm,"eventtime1") <- coxWeibull.lvm(scale=1/100)
    distribution(crm,"eventtime2") <- coxWeibull.lvm(scale=1/100)
    distribution(crm,"censtime") <- coxWeibull.lvm(scale=1/100)
    crm <- eventTime(crm,time~min(eventtime1=1,eventtime2=2,censtime=0),"event")
    crm
}
