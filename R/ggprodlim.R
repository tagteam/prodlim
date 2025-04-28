### ggprodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (14:32) 
## Version: 
## Last-Updated: Apr 28 2025 (13:35) 
##           By: Thomas Alexander Gerds
##     Update #: 110
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' ggplot2::ggplot support for Kaplan-Meier and Aalen-Johansen estimators
##'
##' Important functionality like facet_grid is not yet supported
##' @title ggplot2::ggplot support for Kaplan-Meier and Aalen-Johansen estimators
##' @param x object obtained with \code{\link{prodlim}}.
##' @param xlim Limits for the x-axis.
##' @param ylim Limits for the y-axis.
##' @param x_breaks Breaks for the x-axis.
##' @param y_breaks Breaks for the y-axis.
##' @param position_atrisk Vector of values within xlim. Specifies where numbers at risk should be positioned on the x-axis.
##' @param conf_int Logical. If \code{TRUE} pointwise confidence intervals as a shadow.
##' @param ... passed on to \code{\link{as.data.table.prodlim}}. Can be used to specify 'cause', 'newdata', and 'times'.
##' @return A ggplot2::ggplot object
##' @seealso \code{\link{plot.prodlim}}
##' @examples
##' library(ggplot2)
##' 
##' # Kaplan-Meier and stratified Kaplan-Meier
##' 
##' set.seed(9)
##' ds <- SimSurv(27)
##' 
##' km <- prodlim(Hist(time,event)~1,data = ds)
##' ggprodlim(km)
##' g <- ggprodlim(km)
##' g <- g+geom_step(linewidth=1.5)
##' g + theme(text = element_text(size=20)) + update_geom_defaults("text", list(size=5.5))
##' km1 <- prodlim(Hist(time,event)~X1,data = ds)
##' ggprodlim(km1)
##'
##' ds$group <- factor(sample(1:5,replace=TRUE,size=27),labels=letters[1:5])
##' km2 <- prodlim(Hist(time,event)~group,data = ds)
##' ggprodlim(km2,conf_int=FALSE)
##' 
##' # Aalen-Johansen and stratified Aalen-Johansen
##'
##' set.seed(8)
##' d <- SimCompRisk(27)
##' d$X_group <- factor(sample(1:5,replace=TRUE,size=27),labels=letters[1:5])
##' aj <- prodlim(Hist(time,event)~1,data = d)
##' ggprodlim(aj)
##' ggprodlim(aj,cause=1)
##' ggprodlim(aj,position_atrisk=c(0,5,10))+scale_x_continuous(breaks=c(0,5,10))
##' 
##' ggprodlim(aj)+theme_minimal()+theme(plot.margin=margin(t=0,r=0,b=8,l=0,"line"))
##'
##' # changing colors
##' g+ggplot2::scale_fill_manual(values = 1:2)+ggplot2::scale_color_manual(values=1:2)
##' 
##' aj <- prodlim(Hist(time,event)~X1,data = d)
##' ggprodlim(aj,cause = 1)
##' d$X1 <- factor(d$X1,levels=c("1","0"),labels=c("1","0"))
##' aj <- prodlim(Hist(time,event)~X1,data = d)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
ggprodlim <- function(x,
                      xlim,
                      ylim,
                      y_breaks,
                      x_breaks,
                      position_atrisk,
                      conf_int,
                      ...){
    time <- absolute_risk <- lower <- upper <- cause <- n.risk <- NULL
    # only discrete covariates (for now)
    if (length(x$continuous.predictors)>0) stop("ggprodlim does not deal with continuous.predictors yet. Use plot.prodlim for now.")
    if (length(x$clustervar)>0) stop("ggprodlim does not deal with clustered data yet. Use plot.prodlim for now.")
    covariates <- x$discrete.predictors
    # create data table with plot data
    w <- data.table::as.data.table(x = x,...)
    covariates <- covariates[sapply(covariates,function(c){length(unique(w[[c]]))>1})]
    if (length(covariates) == 0) covariates <- NULL
    if (length(covariates)>0){
        # make sure the covariates are factors
        w[,(covariates) := lapply(.SD,as.factor),.SDcols = covariates]
    }
    outcome_type <- ifelse(match("absolute_risk",names(w),nomatch = 0)>0,"absolute_risk","surv")
    # detect the colour and fill variables
    color_variable <- covariates
    if (length(unique(w$cause))>1) {
        if (length(covariates)>0){
            stop("Cannot plot multiple causes of a stratified Aalen-Johansen estimator in one graph.\nYou should call the function multiple times, first with cause = 1, then with cause =2, etc.")
        }else{
            color_variable <- "cause"
        }
    }
    # create the object
    g <- ggplot2::ggplot(data = data.frame(w),
                         ggplot2::aes(x = time,
                                      y = !! rlang::sym(outcome_type),
                                      fill = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL},
                                      colour = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL}))
    # getting data for numbers at-risk below the graph
    if (missing(position_atrisk)){
        jump_times <- sort(unique(w$time))
        position_atrisk <- seq(min(jump_times),
                               max(jump_times),
                               (max(jump_times)-min(jump_times))/10)
    }
    atrisk_times <- data.table::data.table(time = position_atrisk)
    if (length(covariates)>0){
        atrisk <- unique(w[,c(covariates,"time","n.risk"),with = FALSE])
        strata <- w[,unique(.SD),.SDcols = covariates]
        atrisk_data <- lapply(1:NROW(strata),function(s){
            atrisk_s <- atrisk[strata[s],on = covariates]
            atrisk_s[atrisk_times,on = "time",roll = TRUE]
        })
    }else{
        atrisk <- unique(w[,c("time","n.risk"),with = FALSE])
        atrisk_data <- list(atrisk[atrisk_times,on = "time",roll = TRUE])
    }
    g <- g+ggplot2::geom_step()
    # cannot import pammtools due to a circular imports
    if (match("lower",names(w),nomatch = 0) >0 &&
        ((missing(conf_int) ||
          (length(conf_int)>0 && conf_int != FALSE)))){
        requireNamespace("pammtools")
        g <- g+pammtools::geom_stepribbon(ggplot2::aes(ymin=lower,
                                                       ymax=upper,
                                                       fill = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL},
                                                       colour = !! if (length(color_variable)>0) {rlang::sym(color_variable)}else{NULL}),
                                          linetype = 0,
                                          alpha=0.2,
                                          show.legend = FALSE)
    }
    # axes
    if (missing(ylim)) ylim <- c(0,1)
    if (missing(xlim)) xlim <- c(0,max(w$time))
    ## g <- g+ggplot2::ylim(0,1)
    if (missing(y_breaks)) y_breaks <- seq(ylim[1],ylim[2],abs(ylim[2]-ylim[1])/4)
    ## throws a warning
    g <- g+ggplot2::scale_y_continuous(limits = ylim,breaks = y_breaks,labels = paste0(100*y_breaks,"%"))
    g <- g+ggplot2::coord_cartesian(ylim = ylim,xlim = xlim,clip = 'off')
    ## g <- g+ggplot2::scale_fill_manual(values = grDevices::palette.colors(palette = "Okabe-Ito"))
    ## g <- g+ggplot2::scale_color_manual(values = grDevices::palette.colors(palette = "Okabe-Ito"))
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442")
    g <- g+ggplot2::scale_fill_manual(values = cbbPalette)
    g <- g+ggplot2::scale_color_manual(values = cbbPalette)
    ## g <- g+ggplot2::theme_bw()+ ggplot2::theme(axis.title.x = ggplot2::element_text(vjust=0))
    for (i in 1:length(atrisk_data)){
        atrisk <- atrisk_data[[i]]
        vpos <- 3+(i*2)
        g <- g+ggplot2::geom_text(data = atrisk,
                                  mapping = ggplot2::aes(x = time,
                                                         y = I(0),
                                                         vjust = !!(vpos),
                                                         label = n.risk,
                                                         fill = NULL,
                                                         colour =  !! if (length(covariates)>0){rlang::sym(covariates)}else {NULL}),
                                  show.legend = FALSE)
    }
    # space for atrisk data
    g <- g+ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,length(atrisk_data)+5,1), "lines"))
    g <- g+ggplot2::ylab(ifelse(outcome_type == "absolute_risk","Absolute risk","Survival probability"))+ggplot2::xlab("Time")
    ## g <- g+ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)))
    ## g <- g+ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0)))
    g <- g+ggplot2::annotate("text",x = 0, y = I(0), vjust = ggplot2::unit(3, "lines"), label = "Number at risk")
    g
}
######################################################################
### ggprodlim.R ends here
