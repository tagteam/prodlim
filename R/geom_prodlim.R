geom_prodlim <- function(mapping = NULL, data = NULL,
                         position = "identity", na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE,
                         type = "risk",
                         cause = NULL,
                         conf_int = TRUE,
                         conf_int_alpha = 0.2,
                         percent = TRUE,
                         timeconverter = NULL,
                         times = NULL,
                         cens.code = "0",
                         ...) {
    modify_aes <- function(a, b) {
        if (is.null(a)) return(b)
        la <- as.list(a)
        lb <- as.list(b)
        out <- utils::modifyList(la, lb)
        class(out) <- class(a)
        out
    }
    if (length(cause)>1){
        mapping <- modify_aes(
            mapping,
            ggplot2::aes(
                         group  = ggplot2::after_stat(Event),
                         colour = ggplot2::after_stat(Event),
                         fill   = ggplot2::after_stat(Event)
                     ))
    }
    ggplot2::layer(
                 stat = StatProdlim,
                 geom = GeomProdlim,
                 data = data,
                 mapping = mapping,
                 position = position,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(
                     na.rm = na.rm,
                     type = type,
                     cause = cause,
                     conf_int = conf_int,
                     conf_int_alpha = conf_int_alpha,
                     percent = percent,
                     timeconverter = timeconverter,
                     times = times,
                     cens.code = cens.code,
                     ...
                 )
             )
}

