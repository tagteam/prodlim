## geom_prodlim()
    ## -> defines layer + mapping recipe

## ggplot_build()
    ## -> runs StatProdlim
    ## -> gets returned columns
    ## -> evaluates after_stat(...) using those columns
    ## -> draws geom
#' Kaplan-Meier and Aalen-Johansen curves inside ggplot2
#'
#' Compute and draw Kaplan-Meier or Aalen-Johansen estimators directly
#' within a \code{ggplot2} workflow. The estimator is calculated by
#' calling \code{\link{prodlim}} on the data supplied to the layer.
#'
#' The statistical transformation is performed separately for each
#' \code{ggplot2} group and panel. Consequently, when used together
#' with \code{\link[ggplot2]{facet_grid}} or
#' \code{\link[ggplot2]{facet_wrap}}, the Kaplan-Meier or
#' Aalen-Johansen estimator is computed automatically within each
#' facet.
#'
#' The layer expects the event history to be supplied through the
#' aesthetics \code{x} (time) and \code{event} (event indicator).
#' Competing risks are detected automatically by
#' \code{\link{prodlim}} if the event variable has more than two
#' levels.
#'
#' By default the layer draws a step function representing the
#' Kaplan-Meier survival curve or the cause-specific cumulative
#' incidence function. Optionally, pointwise confidence intervals
#' can be displayed as shaded rectangles following the stepwise
#' structure of the estimator.
#' 
#' `geom_prodlim()` computes Kaplan-Meier or Aalen-Johansen curves inside the
#' statistical transformation. For this reason, grouping variables must be
#' available to the layer **before the stat is evaluated**. This means that
#' in competing risk settings \code{cause} is not available for \code{\link[ggplot2]{facet_grid}}.
#'
#' \preformatted{
#' ggplot(dat, aes(x = time, event = status)) +
#'   geom_prodlim(
#'     mapping = aes(colour = Letters, fill = Letters, group = Letters),
#'     conf_int = TRUE
#'   )
#' }
#'
#' The confidence bands (`geom_rect`) use `fill`, whereas the survival curves
#' (`geom_step`) use `colour`.
#' 
#'
#' @param mapping Set of aesthetic mappings created by
#'   \code{\link[ggplot2]{aes}}. Must include \code{x} (time) and
#'   \code{event}.
#' @param data A data frame containing the variables used in the plot.
#'   If \code{NULL}, the data are inherited from the plot.
#' @param position Position adjustment passed to the underlying
#'   \code{ggplot2} layers.
#' @param na.rm Logical. If \code{TRUE}, missing values are removed
#'   before estimation.
#' @param inherit.aes Logical. If \code{TRUE} (default), the mapping
#'   is combined with the default mapping at the top level of the plot.
#' @param type Character string specifying the type of estimator:
#'   \code{"surv"} for Kaplan-Meier survival probabilities or
#'   \code{"risk"} for cumulative incidence functions.
#' @param cause For competing risks models, specifies which cause
#'   to plot. If \code{NULL}, the first cause is used. If
#'   \code{"stacked"}, a stacked representation of all causes is drawn.
#' @param conf_int Logical. If \code{TRUE}, draw pointwise confidence
#'   intervals.
#' @param conf_int_alpha Transparency level for the confidence interval
#'   shading.
#' @param percent Logical. If \code{TRUE}, probabilities are displayed
#'   on the percentage scale.
#' @param timeconverter Optional character string specifying a time
#'   scale transformation such as \code{"days2years"} or
#'   \code{"months2years"}.
#' @param times Optional vector of time points at which the estimator
#'   should be evaluated.
#' @param ... Additional arguments passed to
#'   \code{\link{prodlim}} or \code{as.data.table.prodlim}.
#'
#' @details
#' The statistical transformation implemented by this layer fits a
#' \code{\link{prodlim}} model of the form
#'
#' \deqn{Hist(time, event) ~ 1}
#'
#' to the data associated with each \code{ggplot2} group. The resulting
#' Kaplan-Meier or Aalen-Johansen estimates are then converted to a
#' step function suitable for plotting with \code{ggplot2}.
#'
#' Confidence intervals are displayed as rectangles spanning the
#' interval between successive event times, mimicking the visual style
#' used by \code{\link{plot.prodlim}}.
#'
#' @return A \code{ggplot2} layer (or list of layers if confidence
#'   intervals are requested).
#'
#' @seealso
#' \code{\link{prodlim}}, \code{\link{plot.prodlim}},
#' \code{\link[ggplot2]{geom_step}}, \code{\link[ggplot2]{facet_grid}}
#'
#' @examples
#' library(prodlim)
#' library(ggplot2)
#' library(riskRegression)
#' library(data.table)
#'
#' ## simulate competing risks data
#' data(Melanoma)
#' ## overall survival
#' data.table::setDT(Melanoma)
#' Melanoma[,combined:=1*(status!=0)]
#' ## estimation within facets
#' 
#' ## Kaplan-Meier curves 
#' ggplot(Melanoma, aes(x = time, event = combined, colour = sex,fill = sex)) +
#'   geom_prodlim(conf_int = TRUE, type="surv")
#' 
#' ## risks instead of survival now with facet_grid
#' ggplot(Melanoma, aes(x = time, event = combined, colour = sex,fill = sex)) +
#'   geom_prodlim(conf_int = TRUE, type="risk")+facet_grid(~epicel)
#'
#' ## Aalen-Johansen curves: cause 1
#' ggplot(dat, aes(x = time, event = event)) + geom_prodlim(conf_int = TRUE,
#'                                             cause="death.malignant.melanoma")
#'
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE) +
#'   facet_grid(1~ Letters)
#'
#' ## Aalen-Johansen curves for selected causes (vector)
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE, cause = 1:2)

#' ## Aalen-Johansen: multiple groups in each panel
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(mapping=aes(fill=after_stat(cause)),
#'   conf_int = TRUE, cause = 1:2)
#'
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE,cause=1:2) +
#'   facet_grid(1~ Letters)
#' 
#' # stacking risks of all causes 
#' ggplot(dat, aes(x = time, event = event))
#' + geom_prodlim(cause="stacked",aes(fill = after_stat(cause)))
#' + facet_grid(1~ Letters)
#'
#' 
#' @export
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
    geom_prodlim_modify_aes <- function(a, b) {
        if (is.null(a)) return(b)
        la <- as.list(a)
        lb <- as.list(b)
        out <- utils::modifyList(la, lb)
        class(out) <- class(a)
        out
    }

    geom_prodlim_drop_aes <- function(mapping, aes_names) {
        if (is.null(mapping)) return(mapping)
        ml <- as.list(mapping)
        for (nm in aes_names) ml[[nm]] <- NULL
        class(ml) <- class(mapping)
        ml
    }

    stacked <- !is.null(cause) && identical(cause, "stacked")
    ## print("mapping")
    ## print(mapping)
    map_list <- if (is.null(mapping)) list() else as.list(mapping)
    has_user_colour <- !is.null(map_list$colour) || !is.null(map_list$color)
    has_user_fill   <- !is.null(map_list$fill)

    col_expr <- NULL
    if (!is.null(map_list$colour)) col_expr <- map_list$colour
    if (is.null(col_expr) && !is.null(map_list$color)) col_expr <- map_list$color

    geom_prodlim_may_inherit <- isTRUE(inherit.aes) && is.null(mapping)

    multi_curve <- (!is.null(cause) &&
                    (identical(cause, "all") ||
                     identical(cause, "stacked") ||
                     (is.numeric(cause) && length(cause) > 1)))

    need_legend <- inherit.aes == TRUE || multi_curve || has_user_colour || has_user_fill || !is.null(map_list$group)

    ## stacked plot: rectangles only
    if (stacked) {
        mapping_rect <- {
            base <- geom_prodlim_modify_aes(
                mapping,
                ggplot2::aes(
                    xmin  = ggplot2::after_stat(xmin),
                    xmax  = ggplot2::after_stat(xmax),
                    ymin  = ggplot2::after_stat(ymin),
                    ymax  = ggplot2::after_stat(ymax),
                    group = interaction(
                        ggplot2::after_stat(geom_prodlim_groupid),
                        ggplot2::after_stat(cause)
                    )
                )
            )

            if (isTRUE(has_user_colour) && is.null(map_list$group) && !is.null(col_expr)) {
                base <- geom_prodlim_modify_aes(
                    base,
                    rlang::inject(ggplot2::aes(group = !!col_expr))
                )
            }

            if (!isTRUE(has_user_fill) && !geom_prodlim_may_inherit) {
                base <- geom_prodlim_modify_aes(
                    base,
                    ggplot2::aes(fill = ggplot2::after_stat(cause))
                )
            }

            base
        }
        return(
            ggplot2::layer(
                stat = StatProdlim,
                data = data,
                mapping = mapping_rect,
                geom = "rect",
                position = position,
                show.legend = need_legend,
                inherit.aes = inherit.aes,
                params = list(
                    na.rm = na.rm,
                    type = type,
                    cause = cause,
                    conf_int = FALSE,
                    percent = percent,
                    timeconverter = timeconverter,
                    times = times,
                    cens.code = cens.code,
                    colour = NA,
                    geom_prodlim_layer = "step",
                    ...
                )
            )
        )
    }

    ## step layer
    mapping_step <- {
        base <- geom_prodlim_drop_aes(mapping, "fill")

        ## if colour defines groups, make that grouping explicit before the stat
        if (isTRUE(has_user_colour) && is.null(map_list$group) && !is.null(col_expr)) {
            base <- geom_prodlim_modify_aes(
                base,
                rlang::inject(ggplot2::aes(group = !!col_expr))
            )
        }

        ## user-supplied colour: keep it unchanged
        if (isTRUE(has_user_colour)) {
            base
        } else {
            ## default for multi-cause output
            if (isFALSE(inherit.aes)){
                if (length(cause)>1) {
                    geom_prodlim_modify_aes(
                        base,
                        ggplot2::aes(group  = interaction(ggplot2::after_stat(geom_prodlim_groupid),ggplot2::after_stat(cause)),
                                     colour  = interaction(ggplot2::after_stat(geom_prodlim_groupid),ggplot2::after_stat(cause)))
                    )
                }else{
                    geom_prodlim_modify_aes(
                        base,
                        ggplot2::aes(group  = ggplot2::after_stat(cause),
                                     colour  = ggplot2::after_stat(cause))
                    )
                }
            }
        }
    }

    layer_step <- ggplot2::layer(
        stat = StatProdlim,
        data = data,
        mapping = mapping_step,
        geom = "step",
        position = position,
        show.legend = need_legend,
        inherit.aes = inherit.aes,
        params = list(
            na.rm = na.rm,
            type = type,
            cause = cause,
            conf_int = FALSE,
            percent = percent,
            timeconverter = timeconverter,
            times = times,
            cens.code = cens.code,
            geom_prodlim_layer = "step",
            ...
        )
    )

    if (!isTRUE(conf_int)) return(layer_step)
    ## CI layer
    mapping_ci <- {
        
        base <- geom_prodlim_modify_aes(
            mapping,
            ggplot2::aes(xmin  = ggplot2::after_stat(xmin),
                         xmax  = ggplot2::after_stat(xmax),
                         ymin  = ggplot2::after_stat(ymin),
                         ymax  = ggplot2::after_stat(ymax)
                         )
        )
        ## print("groups")
        ## print(map_list$group)
        if (isFALSE(inherit.aes)){
            if (length(cause)>1) {
                base <- geom_prodlim_modify_aes(
                    base,
                    ggplot2::aes(group  = interaction(ggplot2::after_stat(geom_prodlim_groupid),ggplot2::after_stat(cause)),
                                 fill  = interaction(ggplot2::after_stat(geom_prodlim_groupid),ggplot2::after_stat(cause)))
                )
            }else{
                base <- geom_prodlim_modify_aes(
                    base,
                    ggplot2::aes(group  = ggplot2::after_stat(cause),
                                 fill  = ggplot2::after_stat(cause))
                )
            }
        }

        ## if colour defines groups, make that grouping explicit before the stat
        if (isTRUE(has_user_colour) && is.null(map_list$group) && !is.null(col_expr)) {
            base <- geom_prodlim_modify_aes(
                base,
                rlang::inject(ggplot2::aes(group = !!col_expr))
            )
        }

        ## user-supplied fill wins; otherwise default fill by cause
        if (!isTRUE(has_user_fill) && !geom_prodlim_may_inherit) {
            base <- geom_prodlim_modify_aes(
                base,
                ggplot2::aes(fill = ggplot2::after_stat(cause))
            )
        }

        base
    }
    ## print("step")
    ## print(mapping_step)
    ## print("ci")
    ## print(mapping_ci)
    layer_ci <- ggplot2::layer(
                             stat = StatProdlim,
                             data = data,
                             mapping = mapping_ci,
                             geom = "rect",
                             position = position,
                             show.legend = need_legend,
                             inherit.aes = inherit.aes,
                             params = list(
                                 na.rm = na.rm,
                                 type = type,
                                 cause = cause,
                                 conf_int = TRUE,
                                 percent = percent,
                                 timeconverter = timeconverter,
                                 times = times,
                                 cens.code = cens.code,
                                 alpha = conf_int_alpha,
                                 colour = NA,
                                 geom_prodlim_layer = "rect",
                                 ...
                             )
                         )

    list(layer_ci, layer_step)
}





