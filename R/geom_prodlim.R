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
#' @param mapping Set of aesthetic mappings created by
#'   \code{\link[ggplot2]{aes}}. Must include \code{x} (time) and
#'   \code{event}.
#' @param data A data frame containing the variables used in the plot.
#'   If \code{NULL}, the data are inherited from the plot.
#' @param position Position adjustment passed to the underlying
#'   \code{ggplot2} layers.
#' @param na.rm Logical. If \code{TRUE}, missing values are removed
#'   before estimation.
#' @param show.legend Logical. Should this layer be included in the
#'   legend?
#' @param inherit.aes Logical. If \code{TRUE} (default), the mapping
#'   is combined with the default mapping at the top level of the plot.
#' @param type Character string specifying the type of estimator:
#'   \code{"surv"} for Kaplan-Meier survival probabilities or
#'   \code{"risk"} for cumulative incidence functions.
#' @param cause For competing risks models, specifies which cause
#'   to plot. If \code{NULL}, the first cause is used. If
#'   \code{"stacked"}, a stacked representation of all causes is drawn.
#' @param legend Legend behaviour: \code{"auto"} hides the legend when only one curve is drawn
#' (Kaplan-Meier or a single cause). Use \code{"show"} or \code{"hide"} to override.
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
#' set.seed(89)
#' dat <- sampleData(89)
#' ## combined endpoint
#' dat[,status:=1*(event!=0)]
#' ## estimation within facets
#' dat[, group:=factor(sample(letters[1:3], .N,replace=TRUE))]
#'
#' ## Kaplan Meier  curves
#' ggplot(dat, aes(x = time, event = status, color=group,group=group)) +
#'   geom_prodlim(conf_int = TRUE)
#'
#' ggplot(dat, aes(x = time, event = status)) +
#'   geom_prodlim(conf_int = TRUE) +
#'   facet_grid(1~ group)
#'
#' ## Aalen-Johansen curves cause 1
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE)
#'
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE) +
#'   facet_grid(1~ group)
#'
#' ## Aalen-Johansen curves both causes
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE,cause=1:2)
#'
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE,cause=1:2) +
#'   facet_grid(1~ group)
#'
#' ggplot(dat, aes(x = time, event = event)) +
#'   geom_prodlim(conf_int = TRUE,cause="stacked") +
#'   facet_grid(1~ group)
#'
#' 
#' @export
geom_prodlim <- function(mapping=NULL, data=NULL,
                         position="identity", na.rm=FALSE, show.legend=NA,
                         inherit.aes=TRUE,
                         type=c("surv","risk"),
                         cause=NULL,
                         legend=c("auto","show","hide"),
                         conf_int=FALSE,
                         conf_int_alpha=0.2,
                         percent=TRUE,
                         timeconverter=NULL,
                         times=NULL,
                         ...) {

    type <- match.arg(type)
    stacked <- (!is.null(cause) && identical(cause, "stacked"))
    legend <- match.arg(legend)
    multi_curve <- (!is.null(cause) && (identical(cause, "all") || identical(cause, "stacked") || (is.numeric(cause) && length(cause) > 1)))
    if (legend == "hide") {
        show_leg_step <- FALSE
        show_leg_ci <- FALSE
    } else if (legend == "show") {
        show_leg_step <- TRUE
        show_leg_ci <- TRUE
    } else {
        ## auto: only show legend when user requested multiple causes
        show_leg_step <- isTRUE(multi_curve)
        show_leg_ci <- isTRUE(multi_curve)
    }


    ##--------------------------------------------------------------
    ## stacked plot: one rect layer (no step curve)
    if (stacked){

        mapping_rect <- .gg_merge_aes(
            mapping,
            ggplot2::aes(
                         xmin = ggplot2::after_stat(xmin),
                         xmax = ggplot2::after_stat(xmax),
                         ymin = ggplot2::after_stat(ymin),
                         ymax = ggplot2::after_stat(ymax),
                         group = ggplot2::after_stat(cause),
                         fill  = ggplot2::after_stat(cause)
                     )
        )

        return(
            ggplot2::layer(
                         stat = StatProdlim,
                         data = data,
                         mapping = mapping_rect,
                         geom = "rect",
                         position = position,
                         show.legend = show_leg_step,
                         inherit.aes = inherit.aes,
                         params = list(
                             na.rm = na.rm,
                             type = type,
                             cause = cause,
                             conf_int = FALSE,
                             percent = percent,
                             timeconverter = timeconverter,
                             times = times,
                             colour = NA,
                             ...
                         )
                     )
        )
    }

    ##--------------------------------------------------------------
    ## step curve layer (filters to .gg_kind == "step")
    mapping_step <- .gg_merge_aes(
        mapping,
        ggplot2::aes(
                     group  = ggplot2::after_stat(cause),
                     colour = ggplot2::after_stat(cause)
                 )
    )

    layer_step <- ggplot2::layer(
                               stat = StatProdlim,
                               data = data,
                               mapping = mapping_step,
                               geom = "step",
                               position = position,
                               show.legend = show_leg_step,
                               inherit.aes = inherit.aes,
                               params = list(
                                   na.rm = na.rm,
                                   type = type,
                                   cause = cause,
                                   conf_int = FALSE,
                                   percent = percent,
                                   timeconverter = timeconverter,
                                   times = times,
                                   .gg_kind_keep = "step",
                                   ...
                               )
                           )

    if (!isTRUE(conf_int))
        return(layer_step)

    ##--------------------------------------------------------------
    ## CI shadow via rects (filters to .gg_kind == "rect")
    mapping_ci <- .gg_merge_aes(
        mapping,
        ggplot2::aes(
                     xmin = ggplot2::after_stat(xmin),
                     xmax = ggplot2::after_stat(xmax),
                     ymin = ggplot2::after_stat(ymin),
                     ymax = ggplot2::after_stat(ymax),
                     group = ggplot2::after_stat(cause),
                     fill  = ggplot2::after_stat(cause)
                 )
    )

    layer_ci <- ggplot2::layer(
                             stat = StatProdlim,
                             data = data,
                             mapping = mapping_ci,
                             geom = "rect",
                             position = position,
                             show.legend = show_leg_ci,
                             inherit.aes = inherit.aes,
                             params = list(
                                 na.rm = na.rm,
                                 type = type,
                                 cause = cause,
                                 conf_int = TRUE,
                                 percent = percent,
                                 timeconverter = timeconverter,
                                 times = times,
                                 alpha = conf_int_alpha,
                                 colour = NA,
                                 .gg_kind_keep = "rect",
                                 ...
                             )
                         )

    list(layer_ci, layer_step)
} 


#' @keywords internal
StatProdlim <- ggplot2::ggproto("StatProdlim", ggplot2::Stat,
                                required_aes = c("x","event"),
                                dropped_aes  = "event",
                                compute_group = function(data, scales,
                                                         type = c("surv","risk"),
                                                         cause = NULL,
                                                         conf_int = FALSE,
                                                         percent = TRUE,
                                                         times = NULL,
                                                         timeconverter = NULL,
                                                         .gg_kind_keep = NULL,
                                                         ...) {
                                    type <- match.arg(type)
                                    if (nrow(data) == 0L) return(data.frame())
                                    if (!"event" %in% names(data))
                                        stop("stat_prodlim: need 'event' variable.")

                                    if (is.logical(data$event)) {
                                        stop(
                                            "geom_prodlim(): `event` must be coded as numeric/integer (0=censor, 1,2,...=events).\n",
                                            "Logical TRUE/FALSE is not supported because ggplot2 treats it as a discrete ",
                                            "aesthetic and may split the data into groups before the statistic is computed.\n\n",
                                            "Use:\n",
                                            "  aes(event = as.integer(status))\n",
                                            "or convert the variable beforehand."
                                        )
                                    }
                                    
                                    dots <- list(...)
                                    if (!is.null(times)) dots$times <- times

                                    .gg_asdt <- function(dots_local){
                                        do.call(data.table::as.data.table,
                                                c(list(fit, surv = (type=="surv"), percent = percent),
                                                  dots_local))
                                    }

                                    stacked  <- (!is.null(cause) && identical(cause, "stacked"))
                                    allcauses <- (!is.null(cause) && identical(cause, "all"))
                                    cause_vec <- (!is.null(cause) && is.numeric(cause) && length(cause) >= 1)

                                    fit <- prodlim::prodlim(prodlim::Hist(x, event) ~ 1, data = data, type = type, ...)
                                    
                                    ## as.data.table: decide what causes to request
                                    ## default for competing risks: cause 1 only (unless user asked for all/stacked)
                                    ## (we only know CR after as.data.table, so we do two-stage logic)
                                    w <- .gg_asdt(dots)

                                    if ("cuminc" %in% names(w) && !"absolute_risk" %in% names(w))
                                        data.table::setnames(w, "cuminc", "absolute_risk")

                                    ## if competing risks and user did not specify cause/all/stacked: default to cause 1
                                    if (("cause" %in% names(w)) && (data.table::uniqueN(w$cause) > 1L) &&
                                        is.null(cause)) {

                                        dots$cause <- 1
                                        w <- .gg_asdt(dots)
                                        if ("cuminc" %in% names(w) && !"absolute_risk" %in% names(w))
                                            data.table::setnames(w, "cuminc", "absolute_risk")
                                    } else {
                                        ## if user requested specific cause(s), apply it/them
                                        if (cause_vec && !allcauses && !stacked) {
                                            dots$cause <- cause
                                            w_try <- try(.gg_asdt(dots), silent = TRUE)
                                            if (inherits(w_try, "try-error")) {
                                                ## fallback: request all causes then subset
                                                dots$cause <- NULL
                                                w <- .gg_asdt(dots)
                                                if ("cause" %in% names(w))
                                                    w <- w[w$cause %in% cause]
                                            } else {
                                                w <- w_try
                                            }
                                            if ("cuminc" %in% names(w) && !"absolute_risk" %in% names(w))
                                                data.table::setnames(w, "cuminc", "absolute_risk")
                                        }
                                    }

                                    ## time conversion
                                    if (!is.null(timeconverter)) {
                                        conv <- c(days2years=1/365.25,
                                                  months2years=1/12,
                                                  days2months=1/30.4368499,
                                                  years2days=365.25,
                                                  years2months=12,
                                                  months2days=30.4368499)
                                        if (!timeconverter %in% names(conv))
                                            stop("stat_prodlim: unknown timeconverter")
                                        w[, time := time * unname(conv[timeconverter])]
                                    }

                                    ## y column
                                    if (type == "surv") {
                                        if ("surv" %in% names(w)) {
                                            w[, .gg_y := surv]
                                        } else if ("absolute_risk" %in% names(w)) {
                                            ## fallback (binary event): surv = 1 - risk
                                            w[, .gg_y := (if (percent) 100 else 1) - absolute_risk]
                                        } else stop("stat_prodlim: cannot find outcome column.")
                                    } else {
                                        ## type == "risk"
                                        if ("absolute_risk" %in% names(w)) {
                                            w[, .gg_y := absolute_risk]
                                        } else if ("surv" %in% names(w)) {
                                            ## fallback (binary event): risk = 1 - surv
                                            w[, .gg_y := (if (percent) 100 else 1) - surv]
                                        } else stop("stat_prodlim: cannot find outcome column.")
                                    }

                                    ## CI columns
                                    if (isTRUE(conf_int) && all(c("lower","upper") %in% names(w))) {
                                        w[, .gg_lower := lower]
                                        w[, .gg_upper := upper]
                                    } else {
                                        w[, .gg_lower := NA_real_]
                                        w[, .gg_upper := NA_real_]
                                    }

                                    ## if we used the fallback transformation between risk and surv, transform CI accordingly
                                    if (isTRUE(conf_int) && all(c("lower","upper") %in% names(w))) {
                                        if (type == "surv" && !("surv" %in% names(w)) && ("absolute_risk" %in% names(w))) {
                                            tmp_lower <- w$.gg_lower
                                            tmp_upper <- w$.gg_upper
                                            w[, .gg_lower := (if (percent) 100 else 1) - tmp_upper]
                                            w[, .gg_upper := (if (percent) 100 else 1) - tmp_lower]
                                        }
                                        if (type == "risk" && !("absolute_risk" %in% names(w)) && ("surv" %in% names(w))) {
                                            tmp_lower <- w$.gg_lower
                                            tmp_upper <- w$.gg_upper
                                            w[, .gg_lower := (if (percent) 100 else 1) - tmp_upper]
                                            w[, .gg_upper := (if (percent) 100 else 1) - tmp_lower]
                                        }
                                    }

                                    ## decide if multiple causes are present in the returned w
                                    multi_cause <- ("cause" %in% names(w)) && (data.table::uniqueN(w$cause) > 1L)

                                    ## make clean step table (do NOT depend on .gg_next)
                                    ## Note: as.data.table.prodlim may return multiple rows at the same time.
                                    ## For right-continuous step functions we keep the LAST row per (cause,time).
                                    if (multi_cause) {
                                        data.table::setorder(w, cause, time)
                                        w_step <- w[, .SD[.N], by=.(cause, time)]
                                        data.table::setorder(w_step, cause, time)
                                        w_step[, cause := factor(cause)]
                                    } else {
                                        data.table::setorder(w, time)
                                        w_step <- w[, .SD[.N], by=time]
                                        data.table::setorder(w_step, time)
                                        ## keep a constant factor for stable grouping/colouring
                                        w_step[, cause := factor(1)]
                                    }
                                    ## keep only the needed columns
                                    w_step <- w_step[, .(time, .gg_y, .gg_lower, .gg_upper, cause)]

                                    ## stacked: return rectangles for stacking (no step curve here)
                                    if (stacked) {
                                        if (!multi_cause) stop("stat_prodlim: cause='stacked' requires multiple causes and absolute_risk.")
                                        if (!("absolute_risk" %in% names(w))) stop("stat_prodlim: cause='stacked' requires absolute_risk.")

                                        dt <- data.table::copy(w)
                                        data.table::setorder(dt, time, cause)
                                        dt[, .gg_ymin := data.table::shift(cumsum(absolute_risk), fill=0), by=time]
                                        dt[, .gg_ymax := cumsum(absolute_risk), by=time]

                                        data.table::setorder(dt, cause, time)
                                        dt[, .gg_next := data.table::shift(time, type="lead"), by=cause]
                                        dt <- dt[!is.na(.gg_next)]
                                        dt[, xmin := time]
                                        dt[, xmax := .gg_next]

                                        return(data.frame(
                                            x = dt$xmin,
                                            xmin = dt$xmin,
                                            xmax = dt$xmax,
                                            ymin = dt$.gg_ymin,
                                            ymax = dt$.gg_ymax,
                                            cause = factor(dt$cause)
                                        ))
                                    }

                                    ## CI rectangles table
                                    w_rect <- data.table::copy(w_step)
                                    if (multi_cause) data.table::setorder(w_rect, cause, time) else data.table::setorder(w_rect, time)

                                    if (multi_cause) {
                                        w_rect[, .gg_next := data.table::shift(time, type="lead"), by=cause]
                                    } else {
                                        w_rect[, .gg_next := data.table::shift(time, type="lead")]
                                    }

                                    w_rect <- w_rect[!is.na(.gg_next)]
                                    w_rect[, xmin := time]
                                    w_rect[, xmax := .gg_next]

                                    ## return both step points and rect intervals, marked by kind
                                    ## ALWAYS return a cause factor (single-curve: cause == 1)
                                    out_step <- data.frame(
                                        x = w_step$time,
                                        y = w_step$.gg_y,
                                        cause = if (multi_cause) factor(w_step$cause) else factor(1),
                                        .gg_kind = "step"
                                    )

                                    out_rect <- data.frame(
                                        xmin = w_rect$xmin,
                                        xmax = w_rect$xmax,
                                        ymin = w_rect$.gg_lower,
                                        ymax = w_rect$.gg_upper,
                                        cause = if (multi_cause) factor(w_rect$cause) else factor(1),
                                        .gg_kind = "rect"
                                    )
                                    ## make columns identical so rbind always works
                                    out_step$xmin <- NA_real_
                                    out_step$xmax <- NA_real_
                                    out_step$ymin <- NA_real_
                                    out_step$ymax <- NA_real_

                                    out_rect$x <- NA_real_
                                    out_rect$y <- NA_real_

                                    out <- rbind(out_step, out_rect)

                                    ## keep only what the geom needs (set by geom_prodlim)
                                    if (!is.null(.gg_kind_keep)) {
                                        out <- out[out$.gg_kind %in% .gg_kind_keep, , drop = FALSE]
                                    }
                                    ## Force stable plotting order for the step layer
                                    out <- out[order(out$cause, out$x, out$.gg_kind), ]
                                    out
                                    out
                                }
                                )

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
                     ...
                 )
             )
}




## merge two aes() objects without losing class
.gg_merge_aes <- function(a, b){
    if (is.null(a)) return(b)
    la <- as.list(a)
    lb <- as.list(b)
    out <- utils::modifyList(la, lb)
    class(out) <- class(a)
    out
}
utils::globalVariables(c(
  "xmin","xmax","ymin","ymax","cause",".gg_kind"
))
