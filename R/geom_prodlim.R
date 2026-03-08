geom_prodlim <- function(mapping = NULL, data = NULL,
                         position = "identity", na.rm = FALSE,
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
    drop_aes <- function(mapping, aes_names) {
        if (is.null(mapping)) return(mapping)
        ml <- as.list(mapping)
        for (nm in aes_names) ml[[nm]] <- NULL
        class(ml) <- class(mapping)
        ml
    }
    # mapping contains aes specified directly inside geom_prodlim(aes(...))
    map_list <- if (is.null(mapping)) list() else as.list(mapping)
    has_user_colour <- !is.null(map_list$colour) || !is.null(map_list$color)
    has_user_fill <- !is.null(map_list$fill)
    has_user_group <- !is.null(map_list$group)
    col_expr <- NULL
    if (!is.null(map_list$colour)) col_expr <- map_list$colour
    if (is.null(col_expr) && !is.null(map_list$color)) col_expr <- map_list$color
    ## need_legend <- isTRUE(multi_cause) || has_user_colour || has_user_fill || has_user_group
    # aes contains x and event which after stat becomes x and y
    mapping_step <- modify_aes(
        ## drop_aes(mapping, "fill"),
        mapping,
        ggplot2::aes(group = ggplot2::after_stat(curve_id),
                     colour = ggplot2::after_stat(curve_id))
    )
    # the Kaplan-Meier or Aalen-Johansen line
    layer_step <- ggplot2::layer(stat = StatProdlim,data = data,mapping = mapping_step,geom = "step",position = position,inherit.aes = inherit.aes,params = list(na.rm = na.rm,type = type,cause = cause,conf_int = FALSE,percent = percent,timeconverter = timeconverter,times = times,cens.code = cens.code,...))
    if (!isTRUE(conf_int)) return(layer_step)
    # after stat aes for geom_rect
    mapping_ci <- modify_aes(mapping,ggplot2::aes(xmin = ggplot2::after_stat(xmin),
                                                  xmax = ggplot2::after_stat(xmax),
                                                  ymin = ggplot2::after_stat(ymin),
                                                  ymax = ggplot2::after_stat(ymax),
                                                  ## fill = ggplot2::after_stat(curve_id),
                                                  ## colour = ggplot2::after_stat(curve_id),
                                                  group = ggplot2::after_stat(curve_id)))
    # confidence shadow layer
    layer_ci <- ggplot2::layer(stat = StatProdlim,
                               data = data,
                               mapping = mapping_ci,
                               geom = "rect",
                               position = position,
                               inherit.aes = inherit.aes,
                               params = list(na.rm = na.rm,type = type,cause = cause,conf_int = TRUE,percent = percent,timeconverter = timeconverter,times = times,cens.code = cens.code,alpha = conf_int_alpha,...))
    list(layer_ci, layer_step)
}
