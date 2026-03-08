#' Product-limit geoms for survival and competing risks curves
#'
#' Draw Kaplan-Meier or Aalen-Johansen curves with optional confidence bands.
#' For competing risks with multiple causes, the default display shows one curve
#' per cause. When `cause = "stacked"`, causes are drawn as stacked filled
#' rectangles instead of step curves with confidence shadows.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data Data frame.
#' @param position Position adjustment.
#' @param na.rm If `FALSE`, missing values are removed with a warning.
#' @param show.legend Logical. Should this layer be included in the legends?
#' @param inherit.aes If `FALSE`, overrides the default aesthetics.
#' @param type Passed to [prodlim::prodlim()], usually `"risk"` or `"surv"`.
#' @param cause Cause(s) for competing risks. Can be a vector of causes or the
#'   string `"stacked"`.
#' @param conf_int Logical. Draw confidence intervals.
#' @param conf_int_alpha Alpha level for confidence shadows.
#' @param percent Logical. Passed to [prodlim::summary.prodlim()].
#' @param timeconverter Optional time conversion string.
#' @param times Optional evaluation times passed to
#'   [prodlim::summary.prodlim()].
#' @param cens.code Censoring code passed to [prodlim::Hist()].
#' @param ... Further arguments passed to [ggplot2::layer()].
#'
#' @details
#' When multiple causes are specified, `fill`/`colour` aesthetics are ignored
#' and replaced by a cause-based mapping so that all causes are shown in a
#' single legend. When `cause = "stacked"`, causes are stacked on top of each
#' other and confidence shadows are not drawn.
#'
#' @examples
#' library(riskRegression)
#' library(data.table)
#' library(ggplot2)
#' data(Melanoma)
#' # Kaplan-Meier
#' ggplot(data = Melanoma,aes(x = time, event = 1*(status != 0)))+geom_prodlim(type = "surv")
#' # stratified Kaplan-Meier inherited aes
#' ggplot(data = Melanoma,aes(x = time, event = 1*(status != 0),fill = sex,color = sex))+geom_prodlim(type = "surv")
#' # stratified Kaplan-Meier geom aes
#' ggplot(data = Melanoma,aes(x = time, event = 1*(status != 0)))+geom_prodlim(aes(fill = sex,color = sex),type = "surv")
#' # facet
#' ggplot(data = Melanoma,aes(x = time, event = 1*(status != 0)))+geom_prodlim(type = "surv")+facet_grid(~sex)
#' # stratified and facet
#' ggplot(data = Melanoma,aes(x = time, event = 1*(status != 0),fill = epicel,color = epicel))+geom_prodlim(type = "surv")+facet_grid(~sex)
#'
#' # Aalen-Johansen
#' ggplot(data = Melanoma,aes(x = time, event = status))+geom_prodlim(type = "surv",cause = 1:2)
#' # stratified Aalen-Johansen inherited aes
#' ggplot(data = Melanoma,aes(x = time, event = status,fill = sex,color = sex))+geom_prodlim(type = "surv")
#' # stratified Aalen-Johansen geom aes
#' ggplot(data = Melanoma,aes(x = time, event = status))+geom_prodlim(aes(fill = sex,color = sex),type = "surv")
#' # facet
#' ggplot(data = Melanoma,aes(x = time, event = status))+geom_prodlim(type = "surv")+facet_grid(~sex)
#' # stratified and facet
#' ggplot(data = Melanoma,aes(x = time, event = status,fill = epicel,color = epicel))+geom_prodlim(type = "surv")+facet_grid(~sex)
#' # stacked
#' ggplot(data = Melanoma,aes(x = time, event = status))+geom_prodlim(cause = "stacked")+facet_grid(~sex)
#' # stacked with cens.code option
#' ggplot(data = Melanoma,aes(x = time, event = event))+geom_prodlim(cens.code="censored",cause = "stacked")+facet_grid(~sex)
#' @export
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

    multi_cause <- identical(cause, "stacked") || length(cause) > 1
    user_mapping <- if (is.null(mapping)) character(0) else names(as.list(mapping))

    if (multi_cause) {
        mapping <- modify_aes(
            mapping,
            ggplot2::aes(
                group  = ggplot2::after_stat(Event),
                colour = ggplot2::after_stat(Event),
                fill   = ggplot2::after_stat(Event)
            )
        )
        if (length(intersect(user_mapping, c("colour", "color", "fill"))) > 0) {
            warning(
                "colour/fill aesthetics are ignored when multiple causes are specified",
                call. = FALSE
            )
        }
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
