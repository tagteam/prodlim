### ggprodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (14:32) 
## Version: 
## Last-Updated: Mar  4 2025 (15:58) 
##           By: Thomas Alexander Gerds
##     Update #: 58
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
##' @param x_at Vector containg the positions for ticks on the x-axis
##' @param pos_atrisk Vector with x positions for numbers at risk 
##' @param ... passed on to \code{\link{as.data.table.prodlim}}. Can be used to specify 'cause', 'newdata', and 'times'.
##' @return A ggplot2::ggplot object
##' @seealso \code{\link{plot.prodlim}}
##' @examples
##' library(ggplot2)
##' d <- SimCompRisk(117)
##' fit1 <- prodlim(Hist(time,event)~1,data = d)
##' ggprodlim(fit1,x_at = c(0,4,8,12))
##' g <- ggprodlim(fit1,x_at = c(0,4,8,12))+ggthemes::theme_solarized()
##' g+theme(plot.margin = margin(t = 0,r = 0,b = 4,l = 0,"line"))
##' g <- g+ggplot2::scale_fill_manual(values = grDevices::palette.colors(palette = "Okabe-Ito"))
##' g <- g+ggplot2::scale_color_manual(values = grDevices::palette.colors(palette = "Okabe-Ito"))
##' d$X1 <- factor(d$X1,levels=c("1","0"),labels=c("1","0"))
##' fit <- prodlim(Hist(time,event)~X1,data = d)
##' ggprodlim(fit,cause = 1,x_at = c(0,6,12))
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
ggprodlim <- function(x,x_at,pos_atrisk,...){
    time <- absolute_risk <- cause <- n.risk <- NULL
    ## cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442")
    if (missing(pos_atrisk)) pos_atrisk <- x_at
    w <- data.table::as.data.table(x = x,...)
    if (length(x$continuous.predictors)>0) stop("ggprodlim does not deal with continuous.predictors yet. Use plot.prodlim until then.")
    covariates <- x$discrete.predictors
    atrisk_times <- data.table::data.table(time = x_at)
    if (length(covariates)>0){
        w[,(covariates) := lapply(.SD,as.factor),.SDcols = covariates]
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
    type <- factor(interaction(length(covariates) > 0,length(unique(w$cause))>1),levels = c("FALSE.FALSE","FALSE.TRUE","TRUE.FALSE","TRUE.TRUE"),
                   labels = c("NoCovaNoCause","NoCovaHasCause","HasCovaNoCause","HasCovaHasCause"))
    # cannot import pammtools due to a circular imports
    requireNamespace("pammtools")
    switch(as.character(type),
           "NoCovaNoCause" = {
               g <- ggplot2::ggplot(data = data.frame(w),ggplot2::aes(x = time,y = absolute_risk))+ggplot2::geom_step()
           },"NoCovaHasCause" = {
               g <- ggplot2::ggplot(data = data.frame(w),ggplot2::aes(x = time,y = absolute_risk,col = cause))
               g <- g+pammtools::geom_stepribbon(ggplot2::aes_string(ymin="lower",ymax="upper",fill = "cause"),linetype = 0,alpha=0.2)
           },"HasCovaNoCause" = {
               g <- ggplot2::ggplot(data = data.frame(w),ggplot2::aes_string(x = "time",y = "absolute_risk",col = covariates,fill = covariates))
               g <- g+pammtools::geom_stepribbon(ggplot2::aes_string(ymin="lower",ymax="upper",fill = covariates),linetype = 0,alpha=0.2)
           },"HasCovaHasCause" = {
               g <- ggplot2::ggplot(data = data.frame(w),ggplot2::aes_string(x = "time",y = "absolute_risk",col = c(covariates,"cause"),fill = c(covariates,"cause")))
               g <- g+pammtools::geom_stepribbon(ggplot2::aes_string(ymin="lower",ymax="upper",fill = c(covariates,"cause")),linetype = 0,alpha=0.2)        
           })
    g <- g+ggplot2::geom_step()
    g <- g+ggplot2::ylim(0,1)
    g <- g+ggplot2::scale_x_continuous(breaks = x_at)
    suppressMessages(g <- g+ggplot2::scale_y_continuous(breaks = seq(0,1,.25),labels = paste0(100*seq(0,1,.25),"%")))
    g <- g+ggplot2::coord_cartesian(ylim = c(0,1),clip = 'off')
    g <- g+ggplot2::scale_fill_manual(values = grDevices::palette.colors(palette = "Okabe-Ito"))
    g <- g+ggplot2::scale_color_manual(values = grDevices::palette.colors(palette = "Okabe-Ito"))
    g <- g+ggplot2::theme_bw()+ ggplot2::theme(axis.title.x = ggplot2::element_text(vjust=0))
    for (i in 1:length(atrisk_data)){
        atrisk <- atrisk_data[[i]]
        vpos <- 7+(i*2)
        if (type %in% c("NoCovaHasCause","NoCovaNoCause")){
            g <- g+ggplot2::geom_text(data = atrisk,
                                      mapping = ggplot2::aes(x = time,y = I(0),vjust = !!(vpos),label = n.risk,fill = NULL,col = NULL),
                                      show.legend = FALSE)
        }else{
            g <- g+ggplot2::geom_text(data = atrisk,
                                      mapping = ggplot2::aes(x = time,y = I(0),vjust = !!(vpos),label = n.risk,fill = NULL),
                                      show.legend = FALSE)
        }
    }
    g <- g+ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,5,1), "lines"))
    g <- g+ggplot2::ylab("Absolute risk")+ggplot2::xlab("Time")
    g <- g+ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)))
    g <- g+ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0)))
    g <- g+ggplot2::annotate("text",x = 0, y = I(0), vjust = 5, label = "Number at risk")
    g
}
######################################################################
### ggprodlim.R ends here
