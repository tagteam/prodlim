### ggprodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (14:32) 
## Version: 
## Last-Updated: mar 10 2026 (14:01) 
##           By: Thomas Alexander Gerds
##     Update #: 428
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' ggplot2 support for Kaplan-Meier and Aalen-Johansen estimators
#'
#' A ggplot2-based alternative to \code{\link{plot.prodlim}} that returns a
#' \code{ggplot} object, enabling ggplot theming, facetting, and composition.
#'
#' @title ggplot2 support for prodlim objects
#' @param x An object obtained with \code{\link{prodlim}}.
#' @param type Either \code{"surv"} or \code{"risk"}. Defaults to \code{x$type}.
#' @param cause For competing risk models: a cause (passed to
#'   \code{\link{as.data.table.prodlim}}), \code{"all"} for all causes, or
#'   \code{"stacked"} for a stacked cumulative incidence plot.
#' @param xlim,ylim Limits for the axes.
#' @param x_breaks,y_breaks Breaks for axes.#
#' @param position_atrisk Numeric vector specifying x-positions where numbers at
#'   risk should be printed. If missing, a default grid is chosen.
#' @param atrisk_title Title shown left of the at-risk table.
#' @param atrisk_labels Optional replacement labels for strata/cause groups in
#'   the at-risk table (in plotting order).
#' @param atrisk_title_line Numeric. Position of the at-risk title below the
#'   plotting region, as an offset from \code{min(ylim)} in units of
#'   \code{diff(ylim)}.
#' @param atrisk_title_step Numeric. Vertical distance from the title to the
#'   first at-risk row, in units of \code{diff(ylim)}.
#' @param atrisk_line_step Numeric. Vertical spacing between at-risk rows, in
#'   units of \code{diff(ylim)}.
#' @param atrisk_bottom_margin Numeric. Bottom plot margin in lines when the
#'   at-risk table is shown.
#' @param facet_panel_spacing_y Optional numeric. Vertical spacing between facet
#'   rows in lines. If \code{NULL} and faceting is used together with an at-risk
#'   table, the default is \code{5 + n_colors}. The same value is used as the
#'   default bottom plot margin below the lowest facet row.
#' @param conf_int Logical. If TRUE and \code{lower}/\code{upper} are available,
#'   draw pointwise confidence intervals.
#' @param conf_int_alpha Alpha for CI shading.
#' @param percent Logical. If TRUE, y axis is shown in percent (0-100).
#' @param timeconverter Optional character, see \code{\link{plot.prodlim}} for
#'   supported values (e.g., "years2months").
#' @param color_var Covariate name must be included in \code{colnames(x$X))} and
#' in competing risk models the word \code{"cause"}.
#' @param facet_formula Formula for facet_wrap. May include covariate names \code{colnames(x$X)} and
#' in competing risk models also the word \code{"cause"}.
#' @param palette Optional vector of colors. If provided, applied via
#'   \code{scale_color_manual}/\code{scale_fill_manual}. Default NULL (do not
#'   impose a palette).
#' @param return_data Logical. If TRUE, return a list with plot and data.
#' @param ... Passed to \code{\link{as.data.table.prodlim}}; commonly \code{newdata},
#'   \code{times}, and for competing risks \code{cause}.
#'
#' @return A ggplot object, or a list if \code{return_data=TRUE}.
#' @seealso \code{\link{plot.prodlim}}, \code{\link{as.data.table.prodlim}}
#'
#' @examples
#' data(Melanoma,package = "riskRegression")
#' # marginal Kaplan-Maier
#' km <- prodlim(Hist(time,status != 0)~1,data = Melanoma)
#' ggprodlim(km)
#' # marginal Aalen-Johansen
#' aj <- prodlim(Hist(time,status)~1,data = Melanoma)
#' ggprodlim(aj)
#' ggprodlim(aj,facet_formula = ~cause)
#'
#' ajsex <- prodlim(Hist(time,status)~sex,data = Melanoma)
#' ggprodlim(ajsex)
#' ggprodlim(ajsex,facet_formula = ~cause)
#' ggprodlim(ajsex,facet_formula = ~sex)
#' ggprodlim(ajsex,color_var = "sex")
#' ggprodlim(ajsex,color_var = "cause")
#'
#' # stacking the causes
#' ggprodlim(ajsex,cause = "stacked")
#' # stacked with two covariates 
#'    ajepisex <- prodlim(Hist(time,status)~sex+epicel,data = Melanoma)
#'    ggprodlim(ajepisex,
#'              cause = "stacked",
#'              facet_formula = sex~epicel,
#'              atrisk_title_line = 0.2)
#'
#' @export
ggprodlim <- function(x,
                      type = x$type,
                      cause,
                      xlim,
                      ylim,
                      y_breaks,
                      x_breaks,
                      position_atrisk,
                      atrisk_title = "Number at risk",
                      atrisk_labels = NULL,
                      atrisk_title_line = 0.14,
                      atrisk_title_step = 0.06,
                      atrisk_line_step = 0.04,
                      atrisk_bottom_margin = NULL,
                      facet_panel_spacing_y = NULL,
                      conf_int = TRUE,
                      conf_int_alpha = 0.2,
                      percent = TRUE,
                      timeconverter,
                      color_var,
                      facet_formula,
                      palette = NULL,
                      return_data = FALSE,
                      ...) {
    stopifnot(inherits(x,"prodlim"))
    # ---- imports & NSE guards
    time <- absolute_risk <- surv <- lower <- upper <- n.risk <- n.lost <- ggprodlim_next_time <- ggprodlim_atrisk_labels <- ggprodlim_curve_id <- ggprodlim_stacked_ymax <- ggprodlim_stacked_ymin <- NULL

    # ---- continuous predictor handling: allow if newdata supplied
    dots <- list(...)
    if ("color_vars" %in% names(dots)){
        stop("ggprodlim: 'color_vars' is not an argument. Change this to 'color_var'.")
    }
    has_newdata <- "newdata" %in% names(dots)
    if (length(x$continuous.predictors) > 0 && !has_newdata) {
        stop("ggprodlim: continuous predictors detected. Please supply 'newdata=' (as in plot.prodlim), or use plot.prodlim().")
    }
    if (has_newdata){
        ggprodlim_newdata <- unique(data.table::as.data.table(dots$newdata))
        dots$newdata <- NULL
        covariates <- data.table::copy(colnames(ggprodlim_newdata))
        has_covariates <- TRUE
    }else{
        covariates <- x$discrete.predictors
        if (length(covariates)>0){
            has_covariates <- TRUE
            ggprodlim_newdata <- data.table::as.data.table(x$X)
        }else{
            has_covariates <- FALSE
        }
    }
    if (!missing(color_var)){
        if (length(color_var) != 1 || !(color_var %in% c(covariates,"cause"))) {
            if (x$model == "competing.risks"){
                stop(paste0("ggprodlim: 'color_var' can be the word 'cause' or must match exactly one of the names of the covariates in the fitted object x$X:\n",
                            paste0(colnames(x$X),collapse = ", ")))
            }else{
                stop(paste0("ggprodlim: 'color_var' must match exactly one of the names covariates in the fitted object x$X:\n",
                            paste0(colnames(x$X),collapse = ", ")))
            }
        }
    }else{
        # may change to 'cause' with multiple causes specified
        color_var <- NULL
    }
    if (!missing(facet_formula)){
        facet_vars <- all.vars(facet_formula)
        if (!all(facet_vars %in% c(covariates,"cause"))){
            stop(paste0("ggprodlim: all variables in facet_formula must be one of the covariates in the fitted object:\n",
                        paste0(colnames(x$X),collapse = ", ")))
        }
    }else{
        facet_vars <- NULL
        facet_formula <- NULL
    }

    # ---- cause handling
    # If user wants "all" causes, we just don't pass 'cause' to as.data.table (so it uses defaults).
    # If user wants "stacked", we need all causes.
    stacked <- FALSE

    if (!missing(cause) && length(cause) > 0) {
        if (identical(cause, "stacked")){
            stacked <- TRUE
            dots$cause <- attr(x$model.response,"states")
        } else {
            if (identical(cause, "all")) {
                dots$cause <- attr(x$model.response,"states")
            }else{
                stopifnot(all(cause %in% attr(x$model.response,"states")))
                dots$cause <- cause
            }
        }
    }
    if (length(dots$cause) == 0){
        requested_causes <- attr(x$model.response,"states")
        cause <- requested_causes
    }else{
        requested_causes <- dots$cause
    }

    # ---- type handling
    if (is.null(type) || length(type) == 0) type <- x$type
    if (!type %in% c("surv","risk")) stop("ggprodlim: 'type' must be 'surv' or 'risk'.")
    if (attr(x$model.response,"model") == "competing.risks"){
        type <- "risk"
    }
    # call summary.prodlim via as.data.table.prodlim

    w <- do.call(data.table::as.data.table,c(list(x = x, surv = (type == "surv"), percent = percent),
                                             dots))

    # normalize column names
    if ("cuminc" %in% names(w) && !"absolute_risk" %in% names(w)) {
        data.table::setnames(w, "cuminc", "absolute_risk")
    }

    # determine outcome column
    outcome_type <- if ("absolute_risk" %in% names(w)) "absolute_risk" else "surv"
    if (type == "risk" && outcome_type == "surv") {
        # as.data.table may still provide surv; convert on the fly if possible
        w[, absolute_risk := 1 - surv]
        outcome_type <- "absolute_risk"
    }

    # check for multiple causes
    has_cause <- ("cause" %in% names(w)) && (length(unique(w$cause)) > 1)
    multi_cause <- has_cause && (length(unique(w$cause))>1)

    if (multi_cause){
        if (!"cause"%in%facet_vars){
            if (length(color_var)>0 && !identical(color_var,"cause")){
                # force cause into facets
                if (length(facet_formula)>0){
                    facet_formula <- reformulate(
                        c(attr(terms(facet_formula), "term.labels"), "cause"),
                        response = if (length(facet_formula) == 3) facet_formula[[2]] else NULL
                    )
                }else{
                    facet_formula <- formula(paste0("~cause"))
                }
                facet_vars <- all.vars(facet_formula)
                ## stop("ggprodlim: when multiple causes are requested 'cause' must appear as color_var or in the facet_formula.")  
                ## facet_vars <- c(facet_vars,"cause")
            }else{
                # if not specified as facet use color to distinguish cause 
                color_var <- "cause"
            }
        }
    }
    if (has_covariates){
        # make sure all fitted covariates appear either in color_var or in facet_vars
        unmapped <- setdiff(names(x$X),c(color_var,facet_vars))
        if (length(unmapped) > 0){
            if (length(color_var) == 0) {
                color_var <- unmapped[[1]]
                if (length(unmapped)>1){
                    stop("ggprodlim: object was fitted with covariates ",paste0(unmapped[[-1]],collapse = ", "),". These must occur either as color_var or in facet_formula.")
                }
            }else{
                if (length(facet_formula)>0){
                    facet_formula <- reformulate(
                        c(attr(terms(facet_formula), "term.labels"), unmapped),
                        response = if (length(facet_formula) == 3) facet_formula[[2]] else NULL
                    )
                }else{
                    facet_formula <- formula(paste0("~",paste(unmapped,collapse = "+")))
                }
                facet_vars <- all.vars(facet_formula)
                ## stop("ggprodlim: object was fitted with covariates ",paste0(unmapped,collapse = ", "),". These must occur either as color_var or in facet_formula.")
            }
        }
    }
    # build a strata column for selection/facetting/atrisk
    # make sure that all covariates are factors
    if (has_covariates){
        w[,(covariates) := lapply(.SD,as.factor),.SDcols = covariates]
    }
    if (multi_cause){
        set(w,j = "cause",value = factor(w[["cause"]],levels = requested_causes))
    }
    if (length(color_var)>0){
        # use newdata to detect the order of the strata
        if (color_var == "cause"){
            levels_color <- levels(w[["cause"]])
        }else{
            if (stacked) {
                stop("ggprodlim: cannot specify color_var for 'stacked' plot.")
            }
            levels_color <- levels(ggprodlim_newdata[[color_var]])
        }
        data.table::setorderv(w,color_var)
        n_colors <- length(levels_color)
        set(w,j = "ggprodlim_curve_id",value = as.integer(w[[color_var]]))
    }else{
        levels_color <- "1"
        n_colors <- 1
        set(w,j = "ggprodlim_curve_id",value = rep(1,NROW(w)))
    }
    
    # ---- default axis limits/breaks
    if (missing(ylim)) {
        ylim <- if (percent) c(0, 100) else c(0, 1)
    }
    if (missing(xlim)) {
        xlim <- c(0, max(w$time, na.rm = TRUE))
    }
    if (missing(y_breaks)) {
        y_breaks <- seq(ylim[1], ylim[2], length.out = 5)
    }
    if (missing(x_breaks)) {
        x_breaks <- ggplot2::waiver()
    }

    # aesthetics
    if (!stacked) {
        if (length(color_var)>0) {
            aes_main <- ggplot2::aes(x = time,
                                     y = !!rlang::sym(outcome_type),
                                     colour = !!rlang::sym(color_var),
                                     fill = !!rlang::sym(color_var),
                                     group = ggprodlim_curve_id)
        }else{
            aes_main <- ggplot2::aes(x = time,
                                     y = !!rlang::sym(outcome_type),
                                     group = ggprodlim_curve_id)
        }
        g <- ggplot2::ggplot(data = w, mapping = aes_main)
    } else {
        g <- NULL
    }

    # ---- confidence intervals: draw as rectangles over intervals
    ggprodlim_add_ci <- function(g, w, xlim, alpha) {
        if (!("lower" %in% names(w) && "upper" %in% names(w))) return(g)
        # prepare interval rectangles per curve
        dt <- data.table::copy(w)
        data.table::setkeyv(dt, c(facet_vars,color_var, "time"))
        dt[, ggprodlim_next_time := data.table::shift(time, type = "lead"), by = setdiff(key(dt),"time")]
        dt <- dt[!is.na(ggprodlim_next_time)]
        dt[, xmin := time]
        dt[, xmax := pmin(ggprodlim_next_time, xlim[2])]
        # clip to xlim
        dt <- dt[xmax > xlim[1]]
        dt[, xmin := pmax(xmin, xlim[1])]
        # map CI
        dt[, ymin := lower]
        dt[, ymax := upper]
        # choose fill aesthetic consistent with main plot
        if (length(color_var)>0) {
            rect_aes <- ggplot2::aes(xmin = xmin,
                                     xmax = xmax,
                                     ymin = ymin,
                                     ymax = ymax,
                                     colour = !!rlang::sym(color_var),
                                     fill = !!rlang::sym(color_var),
                                     group = ggprodlim_curve_id)
        }else{
            rect_aes <- ggplot2::aes(xmin = xmin,
                                     xmax = xmax,
                                     ymin = ymin,
                                     ymax = ymax,
                                     group = ggprodlim_curve_id)
        }
        g + ggplot2::geom_rect(data = dt, mapping = rect_aes,
                               inherit.aes = FALSE, alpha = alpha, colour = NA, show.legend = FALSE)
    }
    if (isTRUE(conf_int)) {
        g <- ggprodlim_add_ci(g, w, xlim = xlim, alpha = conf_int_alpha)
    }
    # set default spacing when facetting is used
    if (length(facet_formula) > 0 && is.null(facet_panel_spacing_y)) {
        facet_panel_spacing_y <- 5 + n_colors
    }
    # ---- main geometry
    if (stacked) {
        if (!has_cause) {
            stop("ggprodlim: cause='stacked' requires a competing risks fit (multiple causes).")
        }
        if (!("absolute_risk" %in% names(w))) {
            stop("ggprodlim: cause='stacked' requires 'absolute_risk' column in data.")
        }

        dt <- data.table::copy(w)

        ## define strata for separate stacked curves/panels, excluding cause
        stack_vars <- intersect(c(covariates, facet_vars), names(dt))
        stack_vars <- setdiff(unique(stack_vars), "cause")

        ## stack vertically at each time within each stratum
        data.table::setorderv(dt, c(stack_vars, "time", "cause"))
        dt[, ggprodlim_stacked_ymax := cumsum(absolute_risk),
           by = c(stack_vars, "time")]
        dt[, ggprodlim_stacked_ymin := data.table::shift(ggprodlim_stacked_ymax, fill = 0),
           by = c(stack_vars, "time")]

        ## create rectangles over time intervals within each cause/stratum
        data.table::setorderv(dt, c(stack_vars, "cause", "time"))
        dt[, ggprodlim_next_time := data.table::shift(time, type = "lead"),
           by = c(stack_vars, "cause")]
        dt <- dt[!is.na(ggprodlim_next_time)]

        dt[, xmin := pmax(time, xlim[1])]
        dt[, xmax := pmin(ggprodlim_next_time, xlim[2])]
        dt <- dt[xmax > xmin]

        g <- ggplot2::ggplot(
                          data = dt,
                          ggplot2::aes(xmin = xmin,
                                       xmax = xmax,
                                       ymin = ggprodlim_stacked_ymin,
                                       ymax = ggprodlim_stacked_ymax,
                                       fill = cause)
                      ) +
            ggplot2::geom_rect(colour = NA)

        if (length(facet_formula)>0){
            g <- g + ggplot2::facet_wrap(facet_formula)
            if (!is.null(facet_panel_spacing_y)) {
                g <- g + ggplot2::theme(panel.spacing.y = grid::unit(facet_panel_spacing_y, "lines"))
            }
        }
    }
    # ---- facetting
    if (!stacked) {
        g <- g + ggplot2::geom_step()
        if (length(facet_formula)>0){
            g <- g + ggplot2::facet_wrap(facet_formula)
            if (!is.null(facet_panel_spacing_y)) {
                g <- g + ggplot2::theme(panel.spacing.y = grid::unit(facet_panel_spacing_y, "lines"))
            }
        }
    }

    # ---- optional palette (do NOT force by default)
    if (!is.null(palette)) {
        g <- g + ggplot2::scale_color_manual(values = palette) +
            ggplot2::scale_fill_manual(values = palette)
    }
    # ---- at-risk table
    atrisk_data <- NULL
    has_atrisk <- FALSE
    if (missing(position_atrisk) || is.numeric(position_atrisk)) {
        has_atrisk <- TRUE
        if (missing(position_atrisk)){
            jump_times <- sort(unique(w$time))
            position_atrisk <- seq(min(jump_times),
                                   max(jump_times),
                                   (max(jump_times)-min(jump_times))/10)
        }else{
            position_atrisk <- unique(position_atrisk)
        }
        atrisk_cols <- intersect(c(facet_vars,color_var, "time", "n.risk","n.lost"),names(w))
        # Numbers at-risk are the same for all causes! 
        if (multi_cause){
            first_cause <- w[1]$cause
            atrisk_data <- w[cause == first_cause, atrisk_cols,with = FALSE]
        }else{
            atrisk_data <- w[, atrisk_cols,with = FALSE]
        }
        if (has_covariates){
            # roll number at risk respecting the covariate strata but ignoring 'cause'
            atrisk_strata_vars <- setdiff(c(facet_vars,color_var),"cause")
            data.table::setorderv(atrisk_data,c(atrisk_strata_vars,"time"))
            # turn into integer for rolling
            atrisk_data[,(atrisk_strata_vars) := lapply(.SD,as.integer),.SDcols = atrisk_strata_vars]
            atrisk_times <- unique(atrisk_data[, atrisk_strata_vars, with = FALSE])
            atrisk_times <- atrisk_times[,data.table::data.table(time = position_atrisk),by = atrisk_strata_vars]
            atrisk_data <- atrisk_data[atrisk_times, on = c(atrisk_strata_vars,"time"), roll = TRUE]
            # now return to factor
            for (asv in atrisk_strata_vars){
                set(atrisk_data,j = asv,value = factor(atrisk_data[[asv]],levels = 1:length(levels(ggprodlim_newdata[[asv]])),labels = levels(ggprodlim_newdata[[asv]])))
            }
        }else{
            # roll number at risk 
            atrisk_times <- data.table::CJ(time = position_atrisk)
            atrisk_data <- atrisk_data[atrisk_times, on = c("time"), roll = TRUE]
        }
    }
    # repeat number at risk for each panel
    if (multi_cause && !identical(color_var,"cause")){
        atrisk_data <- data.table::rbindlist(
                                       lapply(requested_causes, function(value) {
                                           data.table::copy(atrisk_data)[, cause := factor(value,levels = requested_causes)]
                                       })
                                   )
    }
    # Add risk table texts
    if (has_atrisk) {
        # Determine group ordering and labels
        # Build labels: default strata (+ cause if no facet by cause)
        if (length(color_var)>0){
            if (color_var[[1]] == "cause"){
                default_labels <- ""
            }else{
                default_labels <- levels(ggprodlim_newdata[[color_var]])
            }
        }
        if (!is.null(atrisk_labels)) {
            if (length(atrisk_labels) != length(default_labels)) {
                warning("ggprodlim: atrisk_labels length does not match number of curves; ignoring atrisk_labels.")
            } else {
                default_labels <- atrisk_labels
            }
        }
        y_title <- min(ylim) - atrisk_title_line * diff(ylim)
        y_first <- y_title - atrisk_title_step * diff(ylim)
        line_step <- atrisk_line_step * diff(ylim)

        for (i in 1:n_colors){
            if (n_colors>1){
                atrisk_strata <- atrisk_data[atrisk_data[[color_var]] == levels_color[[i]]]
            }else{
                atrisk_strata <- atrisk_data
            }
            set(atrisk_strata,
                j = outcome_type,
                value = y_first - (i - 1) * line_step)

            if (length(color_var) == 0 || identical(color_var,"cause")){
                atrisk_aes <- ggplot2::aes(x = time,
                                           y = !!rlang::sym(outcome_type),
                                           label = n.risk,
                                           color = NULL,
                                           fill = NULL)
            }else{
                atrisk_aes <- ggplot2::aes(x = time,
                                           y = !!rlang::sym(outcome_type),
                                           label = n.risk,
                                           color = !!rlang::sym(color_var),
                                           fill = NULL)
            }
            g <- g + ggplot2::geom_text(data = atrisk_strata,
                                        mapping = atrisk_aes,
                                        inherit.aes = FALSE,
                                        show.legend = FALSE)
        }

        g <- g + ggplot2::annotate("text",
                                   x = xlim[1],
                                   y = y_title,
                                   hjust = 1,
                                   label = atrisk_title)

        if (length(facet_formula) > 0) {
            if (is.null(atrisk_bottom_margin)) {
                atrisk_bottom_margin <- facet_panel_spacing_y
            }
        } else {
            if (is.null(atrisk_bottom_margin)) {
                if (identical(color_var,"cause")){
                    atrisk_bottom_margin <- 9
                }else{
                    atrisk_bottom_margin <- n_colors * 2 + 5
                }
            }
        }
        g <- g + ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, atrisk_bottom_margin, 1), "lines"))
        # space for atrisk data
        g <- g + ggplot2::coord_cartesian(ylim = ylim, xlim = xlim, clip = "off")
        # ---- scales
        if (isTRUE(percent)) {
            g <- g + ggplot2::scale_y_continuous(breaks = y_breaks, labels = function(z) paste0(z, "%"))
        } else {
            g <- g + ggplot2::scale_y_continuous(breaks = y_breaks)
        }
        g <- g + ggplot2::scale_x_continuous(breaks = x_breaks)
    }else{
        # ---- scales
        if (isTRUE(percent)) {
            g <- g + ggplot2::scale_y_continuous(limits = ylim, breaks = y_breaks, labels = function(z) paste0(z, "%"))
        } else {
            g <- g + ggplot2::scale_y_continuous(limits = ylim, breaks = y_breaks)
        }
        g <- g + ggplot2::scale_x_continuous(limits = xlim, breaks = x_breaks)
    }
    # ---- labels
    ylab <- if (outcome_type == "absolute_risk") "Absolute risk" else "Survival probability"
    g <- g + ggplot2::xlab("Time")
    g <- g + ggplot2::ylab(ylab)
    # ---- legend title
    if (isTRUE(return_data)) {
        return(list(plot = g, data = w, atrisk = atrisk_data))
    }
    g
}


######################################################################
### ggprodlim.R ends here
