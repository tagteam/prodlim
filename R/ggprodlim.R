### ggprodlim.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2025 (14:32) 
## Version: 
## Last-Updated: mar  8 2026 (06:57) 
##           By: Thomas Alexander Gerds
##     Update #: 133
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
#'   For competing risks, \code{"risk"} is appropriate.
#' @param cause For competing risk models: a cause (passed to
#'   \code{\link{as.data.table.prodlim}}), \code{"all"} for all causes, or
#'   \code{"stacked"} for a stacked cumulative incidence plot.
#' @param select Optional integer vector selecting which curves to keep, after
#'   curves are identified by unique strata (and cause, if present).
#' @param xlim,ylim Limits for the axes.
#' @param x_breaks,y_breaks Breaks for axes.
#' @param position_atrisk Numeric vector specifying x-positions where numbers at
#'   risk should be printed. If missing, a default grid is chosen.
#' @param atrisk_title Title shown left of the at-risk table.
#' @param atrisk_labels Optional replacement labels for strata/cause groups in
#'   the at-risk table (in plotting order).
#' @param atrisk_show_censored Logical. If TRUE and the data provides
#'   \code{n.cens}, show it in parentheses next to \code{n.risk}.
#' @param atrisk_min Integer. If >0, drop curve segments after n.risk < atrisk_min.
#' @param conf_int Logical. If TRUE and \code{lower}/\code{upper} are available,
#'   draw pointwise confidence intervals.
#' @param conf_int_alpha Alpha for CI shading.
#' @param percent Logical. If TRUE, y axis is shown in percent (0-100).
#' @param timeconverter Optional character, see \code{\link{plot.prodlim}} for
#'   supported values (e.g., "years2months").
#' @param facet One of \code{"none"}, \code{"cause"}, \code{"strata"},
#'   \code{"cause_x_strata"}. When multiple causes and multiple strata exist,
#'   facetting avoids the "one graph only" limitation of the old implementation.
#' @param palette Optional vector of colors. If provided, applied via
#'   \code{scale_color_manual}/\code{scale_fill_manual}. Default NULL (do not
#'   impose a palette).
#' @param return_data Logical. If TRUE, return a list with plot and data.
#' @param ... Passed to \code{\link{as.data.table.prodlim}}; commonly \code{newdata},
#'   \code{times}, and for competing risks \code{cause}.
#'
#' @return A ggplot object, or a list if \code{return_data=TRUE}.
#' @seealso \code{\link{plot.prodlim}}, \code{\link{as.data.table.prodlim}}
#' @export
ggprodlim <- function(x,
                      type = x$type,
                      cause,
                      select = NULL,
                      xlim,
                      ylim,
                      y_breaks,
                      x_breaks,
                      position_atrisk,
                      atrisk_title = "Number at risk",
                      atrisk_labels = NULL,
                      atrisk_show_censored = FALSE,
                      atrisk_min = 0L,
                      conf_int = TRUE,
                      conf_int_alpha = 0.2,
                      percent = TRUE,
                      timeconverter,
                      facet = c("none","cause","strata","cause_x_strata"),
                      palette = NULL,
                      return_data = FALSE,
                      ...) {

  # ---- imports & NSE guards
  time <- absolute_risk <- surv <- lower <- upper <- cause <- n.risk <- n.cens <- NULL
  facet <- match.arg(facet)

  # ---- continuous predictor handling: allow if newdata supplied
  dots <- list(...)
  has_newdata <- "newdata" %in% names(dots)

  if (length(x$continuous.predictors) > 0 && !has_newdata) {
    stop("ggprodlim: continuous predictors detected. Please supply 'newdata=' (as in plot.prodlim), or use plot.prodlim().")
  }

  # ---- determine time conversion factor (if any)
  .time_factor <- 1
  .time_label  <- NULL
  if (!missing(timeconverter) && length(timeconverter) > 0 && !is.null(timeconverter)) {
    conv <- c(
      days2years   = 1/365.25,
      months2years = 1/12,
      days2months  = 1/30.4368499,
      years2days   = 365.25,
      years2months = 12,
      months2days  = 30.4368499
    )
    if (!timeconverter %in% names(conv)) {
      stop("ggprodlim: unknown timeconverter. Supported: ",
           paste(names(conv), collapse = ", "))
    }
    .time_factor <- unname(conv[timeconverter])
    .time_label  <- timeconverter
  }

  # ---- cause handling
  # If user wants "all" causes, we just don't pass 'cause' to as.data.table (so it uses defaults).
  # If user wants "stacked", we need all causes.
  stacked <- FALSE
  if (!missing(cause) && length(cause) > 0) {
    if (identical(cause, "stacked")) {
      stacked <- TRUE
      # do NOT pass cause to as.data.table, we need all
    } else if (!identical(cause, "all")) {
      dots$cause <- cause
    }
  }

  # ---- type handling
  if (is.null(type) || length(type) == 0) type <- x$type
  if (!type %in% c("surv","risk")) stop("ggprodlim: 'type' must be 'surv' or 'risk'.")

  # ---- collect plot data
  w <- data.table::as.data.table(x = x, surv = (type == "surv"), percent = percent, dots)

  # normalize column names
  if ("cuminc" %in% names(w) && !"absolute_risk" %in% names(w)) {
    data.table::setnames(w, "cuminc", "absolute_risk")
  }

  # apply time conversion
  if (.time_factor != 1) {
    w[, time := time * .time_factor]
    if (!missing(xlim)) xlim <- xlim * .time_factor
    if (!missing(x_breaks)) x_breaks <- x_breaks * .time_factor
    if (!missing(position_atrisk)) position_atrisk <- position_atrisk * .time_factor
  }

  # determine outcome column
  outcome_type <- if ("absolute_risk" %in% names(w)) "absolute_risk" else "surv"
  if (type == "risk" && outcome_type == "surv") {
    # as.data.table may still provide surv; convert on the fly if possible
    w[, absolute_risk := 1 - surv]
    outcome_type <- "absolute_risk"
  }

  # ---- identify stratification columns
  covariates <- x$discrete.predictors
  if (length(covariates) > 0) {
    covariates <- covariates[covariates %in% names(w)]
    if (length(covariates) > 0) {
      covariates <- covariates[sapply(covariates, function(cn) length(unique(w[[cn]])) > 1)]
    }
  }
  if (length(covariates) == 0) covariates <- NULL

  # ensure covariates are factors for stable ordering/legend labels
  if (!is.null(covariates)) {
    w[, (covariates) := lapply(.SD, as.factor), .SDcols = covariates]
  }

  has_cause <- ("cause" %in% names(w)) && (length(unique(w$cause)) > 1)

    # ---- build a unified group id for selection/facetting/atrisk
    # strata label (for legend and atrisk)
    if (!is.null(covariates)) {
        w[, .strata_label := do.call(paste, c(.SD, sep = ", ")), .SDcols = covariates]
    } else {
        w[, .strata_label := ""]
    }
    if (has_cause) {
        w[, .cause_label := as.character(cause)]
    } else {
        w[, .cause_label := "1"]
    }
    w[, .curve_id := interaction(.strata_label, .cause_label, drop = TRUE, lex.order = TRUE)]

  # selection
  if (!is.null(select)) {
    u <- levels(w$.curve_id)
    keep <- u[select]
    keep <- keep[!is.na(keep)]
    w <- w[.curve_id %in% keep]
  }

  # atrisk_min: truncate curves
  if (!is.null(atrisk_min) && atrisk_min > 0 && "n.risk" %in% names(w)) {
    w <- w[n.risk >= atrisk_min]
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
    x_breaks <- waiver()
  }

  # ---- base ggplot mapping
  # Color/fill mapping strategy:
  # - If stacked: fill = cause, and no color.
  # - Else: color by strata (if multiple), and optionally facet by cause.
  #
  # facet default:
  # - if user requested explicitly, use it
  # - else: if has_cause && covariates not NULL -> facet by cause
  # - if has_cause && covariates NULL -> no facet, color by cause
  # - else: no facet

  if (facet == "none") {
    if (has_cause && !is.null(covariates)) {
      # sensible default for multi-cause + strata
      facet_eff <- "cause"
    } else {
      facet_eff <- "none"
    }
  } else {
    facet_eff <- facet
  }

  # aesthetics
  if (stacked) {
    aes_main <- ggplot2::aes(x = time, ymin = .ymin, ymax = .ymax, fill = .cause_label, group = .curve_id)
  } else if (has_cause && is.null(covariates)) {
    aes_main <- ggplot2::aes(x = time, y = !!rlang::sym(outcome_type), colour = .cause_label, fill = .cause_label, group = .curve_id)
  } else if (!is.null(covariates)) {
    aes_main <- ggplot2::aes(x = time, y = !!rlang::sym(outcome_type), colour = .strata_label, fill = .strata_label, group = .curve_id)
  } else {
    aes_main <- ggplot2::aes(x = time, y = !!rlang::sym(outcome_type), group = .curve_id)
  }

  g <- ggplot2::ggplot(data = w, mapping = aes_main)

  # ---- confidence intervals (no pammtools): draw as rectangles over intervals
  .add_ci_rect <- function(g, w, xlim, alpha) {
    if (!("lower" %in% names(w) && "upper" %in% names(w))) return(g)

    # prepare interval rectangles per curve
    dt <- data.table::as.data.table(w)
    data.table::setorder(dt, .curve_id, time)
    dt[, .t_next := data.table::shift(time, type = "lead"), by = .curve_id]
    dt <- dt[!is.na(.t_next)]
    dt[, xmin := time]
    dt[, xmax := pmin(.t_next, xlim[2])]
    # clip to xlim
    dt <- dt[xmax > xlim[1]]
    dt[, xmin := pmax(xmin, xlim[1])]

    # map CI in percent if needed
    if (percent) {
      dt[, ymin := lower]
      dt[, ymax := upper]
    } else {
      dt[, ymin := lower]
      dt[, ymax := upper]
    }

    # choose fill aesthetic consistent with main plot
    if (has_cause && is.null(covariates)) {
      rect_aes <- ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = .cause_label, group = .curve_id)
    } else if (!is.null(covariates)) {
      rect_aes <- ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = .strata_label, group = .curve_id)
    } else {
      rect_aes <- ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = .curve_id)
    }

    g + ggplot2::geom_rect(data = dt, mapping = rect_aes,
                           inherit.aes = FALSE, alpha = alpha, colour = NA, show.legend = FALSE)
  }

  if (isTRUE(conf_int)) {
    g <- .add_ci_rect(g, w, xlim = xlim, alpha = conf_int_alpha)
  }

  # ---- main geometry
  if (stacked) {
    # compute stacked ymin/ymax per time & strata
    if (!has_cause) {
      stop("ggprodlim: cause='stacked' requires a competing risks fit (multiple causes).")
    }
    if (!("absolute_risk" %in% names(w))) {
      stop("ggprodlim: cause='stacked' requires 'absolute_risk' column in data.")
    }
    dt <- data.table::as.data.table(w)
    data.table::setorder(dt, .strata_label, time, .cause_label)

    # cumulative sum by time within strata
    dt[, .ymin := data.table::shift(cumsum(absolute_risk), fill = 0), by = .( .strata_label, time )]
    dt[, .ymax := cumsum(absolute_risk), by = .( .strata_label, time )]

    # convert to rectangles over time intervals, per cause (stepwise stacked)
    data.table::setorder(dt, .strata_label, .cause_label, time)
    dt[, .t_next := data.table::shift(time, type = "lead"), by = .( .strata_label, .cause_label )]
    dt <- dt[!is.na(.t_next)]
    dt[, xmin := time]
    dt[, xmax := pmin(.t_next, xlim[2])]
    dt <- dt[xmax > xlim[1]]
    dt[, xmin := pmax(xmin, xlim[1])]

    g <- ggplot2::ggplot(data = dt,
                         ggplot2::aes(xmin = xmin, xmax = xmax, ymin = .ymin, ymax = .ymax,
                                     fill = .cause_label, group = interaction(.strata_label, .cause_label, drop = TRUE))) +
      ggplot2::geom_rect(colour = NA)

    # facet by strata if multiple strata
    if (!is.null(covariates) && length(unique(dt$.strata_label)) > 1) {
      g <- g + ggplot2::facet_wrap(~ .strata_label)
    }

  } else {
    g <- g + ggplot2::geom_step()
  }

  # ---- facetting
  if (!stacked) {
    if (facet_eff == "cause" && has_cause) {
      g <- g + ggplot2::facet_wrap(~ .cause_label)
    } else if (facet_eff == "strata" && !is.null(covariates) && length(unique(w$.strata_label)) > 1) {
      g <- g + ggplot2::facet_wrap(~ .strata_label)
    } else if (facet_eff == "cause_x_strata" && has_cause && !is.null(covariates)) {
      g <- g + ggplot2::facet_grid(.cause_label ~ .strata_label)
    }
  }

  # ---- scales
  if (isTRUE(percent)) {
    g <- g + ggplot2::scale_y_continuous(limits = ylim, breaks = y_breaks, labels = function(z) paste0(z, "%"))
  } else {
    g <- g + ggplot2::scale_y_continuous(limits = ylim, breaks = y_breaks)
  }
  g <- g + ggplot2::scale_x_continuous(limits = xlim, breaks = x_breaks)

  # ---- optional palette (do NOT force by default)
  if (!is.null(palette)) {
    g <- g + ggplot2::scale_color_manual(values = palette) +
      ggplot2::scale_fill_manual(values = palette)
  }

  # ---- at-risk table
  atrisk_data <- NULL
  if (!missing(position_atrisk) && length(position_atrisk) >= 1) {
    if (!("n.risk" %in% names(w))) {
      warning("ggprodlim: no 'n.risk' column available; skipping at-risk table.")
    } else {
      atrisk_times <- data.table::data.table(time = position_atrisk)
      atrisk_cols <- c(".curve_id", ".strata_label", ".cause_label", "time", "n.risk")
      if ("n.cens" %in% names(w)) atrisk_cols <- c(atrisk_cols, "n.cens")

      atrisk <- unique(data.table::as.data.table(w)[, ..atrisk_cols])
      data.table::setorder(atrisk, .curve_id, time)

      atrisk_data <- lapply(split(atrisk, atrisk$.curve_id), function(dti) {
        dti[atrisk_times, on = "time", roll = TRUE]
      })
    }
  } else if (missing(position_atrisk)) {
    # default grid (keeps old behavior)
    jump_times <- sort(unique(w$time))
    if (length(jump_times) > 0) {
      position_atrisk <- seq(min(jump_times), max(jump_times), length.out = 6)
      position_atrisk <- unique(position_atrisk)
      atrisk_times <- data.table::data.table(time = position_atrisk)
      if ("n.risk" %in% names(w)) {
        atrisk_cols <- c(".curve_id", ".strata_label", ".cause_label", "time", "n.risk")
        if ("n.cens" %in% names(w)) atrisk_cols <- c(atrisk_cols, "n.cens")
        atrisk <- unique(data.table::as.data.table(w)[, ..atrisk_cols])
        data.table::setorder(atrisk, .curve_id, time)
        atrisk_data <- lapply(split(atrisk, atrisk$.curve_id), function(dti) {
          dti[atrisk_times, on = "time", roll = TRUE]
        })
      }
    }
  }

  # Add risk table texts
  if (!is.null(atrisk_data) && length(atrisk_data) > 0) {

    # Determine group ordering and labels
    curve_ids <- names(atrisk_data)
    if (is.null(curve_ids)) curve_ids <- levels(w$.curve_id)

    # Build labels: default strata (+ cause if no facet by cause)
    default_labels <- vapply(atrisk_data, function(dt) {
      lab <- unique(dt$.strata_label)
      if (length(lab) == 0) lab <- ""
      if (has_cause && facet_eff == "none" && !is.null(unique(dt$.cause_label))) {
        cl <- unique(dt$.cause_label)
        paste0(lab, " / cause=", cl)
      } else {
        lab
      }
    }, character(1))

    if (!is.null(atrisk_labels)) {
      if (length(atrisk_labels) != length(default_labels)) {
        warning("ggprodlim: atrisk_labels length does not match number of curves; ignoring atrisk_labels.")
      } else {
        default_labels <- atrisk_labels
      }
    }

    # Compose n.risk (and censored) label per time point
    for (i in seq_along(atrisk_data)) {
      dt <- atrisk_data[[i]]
      # text label column
      if (atrisk_show_censored && "n.cens" %in% names(dt)) {
        dt[, .label := paste0(n.risk, " (", n.cens, ")")]
      } else {
        dt[, .label := as.character(n.risk)]
      }
      vpos <- 3 + (i * 1.7)
      g <- g + ggplot2::geom_text(
        data = dt,
        mapping = ggplot2::aes(x = time, y = I(ylim[1]), vjust = vpos, label = .label),
        inherit.aes = FALSE,
        show.legend = FALSE
      )
      # left label
      g <- g + ggplot2::annotate(
        "text",
        x = xlim[1],
        y = I(ylim[1]),
        hjust = 0,
        vjust = vpos,
        label = default_labels[i]
      )
    }

    # title line
    g <- g + ggplot2::annotate(
      "text",
      x = xlim[1],
      y = I(ylim[1]),
      hjust = 0,
      vjust = 3,
      label = atrisk_title
    )

    # space for atrisk data
    g <- g + ggplot2::coord_cartesian(ylim = ylim, xlim = xlim, clip = "off")
    g <- g + ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, length(atrisk_data) * 1.7 + 5, 1), "lines"))
  }

  # ---- labels
  ylab <- if (outcome_type == "absolute_risk") "Absolute risk" else "Survival probability"
  if (!is.null(.time_label)) {
    g <- g + ggplot2::xlab(paste0("Time (", .time_label, ")"))
  } else {
    g <- g + ggplot2::xlab("Time")
  }
  g <- g + ggplot2::ylab(ylab)

  if (isTRUE(return_data)) {
    return(list(plot = g, data = w, atrisk = atrisk_data))
  }
  g
}


######################################################################
### ggprodlim.R ends here
