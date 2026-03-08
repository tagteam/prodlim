#' @keywords internal
stat_prodlim <- function(mapping=NULL, data=NULL, geom="step",
                         position="identity", na.rm=FALSE, show.legend=NA,
                         inherit.aes=TRUE,
                         type=c("surv","risk"),
                         cause=NULL,
                         legend=c("auto","show","hide"),
                         conf_int=FALSE,
                         percent=TRUE,
                         timeconverter=NULL,
                         times=NULL,
                         cens.code,
                         ...) {
    ggplot2::layer(
                 stat = StatProdlim,
                 data = data,
                 mapping = mapping,
                 geom = geom,
                 position = position,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(
                     na.rm = na.rm,
                     type = match.arg(type),
                     cause = cause,
                     conf_int = conf_int,
                     percent = percent,
                     timeconverter = timeconverter,
                     times = times,
                     cens.code = cens.code,
                     ...)
             )
}
