stat_prodlim_setup_data <- function(data, params) {
    data <- data.table::as.data.table(data)
    nuniq <- function(x) data.table::uniqueN(x[!is.na(x)])
    browser(skipCalls=1L)
    # 1) normalize colour/fill so inherited aes also work
    if ("colour" %in% names(data) && !"fill" %in% names(data)) {
        data$fill <- data$colour
    }
    if ("fill" %in% names(data) && !"colour" %in% names(data)) {
        data$colour <- data$fill
    }
    if ("fill" %in% names(data) && "colour" %in% names(data)) {
        miss_fill <- is.na(data$fill)
        if (any(miss_fill)) data$fill[miss_fill] <- data$colour[miss_fill]
        miss_col <- is.na(data$colour)
        if (any(miss_col)) data$colour[miss_col] <- data$fill[miss_col]
    }
    # 2) determine strata from actual evaluated columns
    #    precedence: colour -> fill -> linetype -> shape -> group
    if (!"prodlim_strata" %in% names(data)) {
        strata_source <- NULL
        for (nm in c("colour", "fill", "linetype", "shape")) {
            if (nm %in% names(data) && nuniq(data[[nm]]) > 1L) {
                strata_source <- nm
                break
            }
        }
        if (is.null(strata_source) && "group" %in% names(data) && nuniq(data$group) > 1L) {
            strata_source <- "group"
        }
        if (!is.null(strata_source)) {
            data$prodlim_strata <- data[[strata_source]]
        }
    }
    data
}
StatProdlim <- ggplot2::ggproto(
                            "StatProdlim",
                            ggplot2::Stat,
                            required_aes = c("x", "event"),
                            default_aes = ggplot2::aes(prodlim_strata = after_stat(curve_id)),
                            extra_params = c(
                                "na.rm", "type", "cause", "conf_int", "percent",
                                "times", "timeconverter", "cens.code", "conf_int_alpha"
                            ),
                            setup_data = stat_prodlim_setup_data,
                            compute_panel = function(data, scales,
                                                     type = "risk",
                                                     cause = NULL,
                                                     conf_int = FALSE,
                                                     percent = TRUE,
                                                     times = NULL,
                                                     timeconverter = NULL,
                                                     cens.code = "0",
                                                     conf_int_alpha = 0.2) {

    if (nrow(data) == 0L) return(data.frame())

    if (!"event" %in% names(data)) {
      stop("stat_prodlim: need 'event' variable.")
    }
    if (is.logical(data$event)) {
      stop("geom_prodlim(): `event` must be coded as numeric/integer or factor/character, not logical.")
    }

    data <- data.table::data.table(data)

    ignore_formula <- c(
      "x", "y", "xmin", "xmax", "ymin", "ymax",
      "PANEL", "group", "event", "alpha", "size",
      "linewidth", "linetype", "colour", "color", "fill",
      "component", "curve_id", "prodlim_strata"
    )

    model_cols <- setdiff(names(data), ignore_formula)

    strata_col <- NULL
    if ("prodlim_strata" %in% names(data)) {
      vals <- unique(data$prodlim_strata)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0L) strata_col <- "prodlim_strata"
    }

    rhs_cols <- c(model_cols, strata_col)
    rhs_cols <- unique(rhs_cols[nzchar(rhs_cols)])

    ff <- stats::formula(
      paste0(
        "prodlim::Hist(x, event, cens.code = '", cens.code, "') ~ ",
        if (length(rhs_cols) == 0) "1" else paste(rhs_cols, collapse = "+")
      )
    )

    fit <- prodlim::prodlim(formula = ff, data = data, type = type)

    summary_args <- list(
      x = fit,
      surv = identical(type, "surv"),
      percent = percent
    )

    multi_cause <- FALSE
    if (attr(fit$model.response, "model") == "competing.risks") {
      if (length(cause) == 0) {
        requested_causes <- attr(fit$model.response, "states")[[1]]
      } else {
        cause <- intersect(cause, attr(fit$model.response, "states"))
        if (length(cause) == 0) {
          stop(
            paste0(
              "Cannot find cause(s) ",
              paste0(cause, collapse = ", "),
              ". The values of data$event are ",
              paste0(attr(fit$model.response, "states"), collapse = ", ")
            )
          )
        }
        requested_causes <- cause
      }
      summary_args$cause <- requested_causes
    } else {
      requested_causes <- 1
    }

    if (length(times) > 0) summary_args$times <- times

    w <- do.call(data.table::as.data.table, summary_args)

    if (!("cause" %in% names(w))) {
      data.table::set(w, j = "cause", value = factor(rep(1, nrow(w))))
    } else {
      multi_cause <- data.table::uniqueN(w$cause) > 1L
      if (!multi_cause) {
        w[, cause := factor(requested_causes)]
      } else {
        order_cols <- c(if (length(rhs_cols)) rhs_cols else NULL, "cause", "time")
        data.table::setorderv(w, order_cols)
        w[, cause := factor(cause, levels = requested_causes)]
      }
    }

    if ("cuminc" %in% names(w) && !"absolute_risk" %in% names(w)) {
      data.table::setnames(w, "cuminc", "absolute_risk")
    }

    if (!is.null(timeconverter)) {
      conv <- c(
        days2years = 1 / 365.25,
        months2years = 1 / 12,
        days2months = 1 / 30.4368499,
        years2days = 365.25,
        years2months = 12,
        months2days = 30.4368499
      )
      if (!timeconverter %in% names(conv)) {
        stop("stat_prodlim: unknown timeconverter")
      }
      data.table::set(w, j = "time", value = w$time * unname(conv[timeconverter]))
    }

    y_name <- intersect(c("surv", "absolute_risk"), names(w))

    curve_group_cols <- c()
    if (!is.null(strata_col) && strata_col %in% names(w)) {
      curve_group_cols <- c(curve_group_cols, strata_col)
    }
    if (isTRUE(multi_cause)) {
      curve_group_cols <- c(curve_group_cols, "cause")
    }

    if (length(curve_group_cols) == 0) {
      data.table::setorderv(w, "time")
      data.table::set(w, j = "curve_id", value = factor(rep("1", nrow(w))))
    } else {
      data.table::setorderv(w, c(curve_group_cols, "time"))
      curve_labels <- as.character(w[[curve_group_cols[[1]]]])
      if (length(curve_group_cols) > 1) {
        for (j in 2:length(curve_group_cols)) {
          curve_labels <- paste(curve_labels, w[[curve_group_cols[[j]]]], sep = ", ")
        }
      }
      data.table::set(w, j = "curve_id", value = factor(curve_labels))
    }

    step_dt <- data.table::copy(w)
    data.table::setnames(step_dt, y_name, "y")
    data.table::setnames(step_dt, "time", "x")
    step_dt[, component := "step"]
    step_dt[, xmin := NA_real_]
    step_dt[, xmax := NA_real_]
    step_dt[, ymin := NA_real_]
    step_dt[, ymax := NA_real_]
    step_dt[, group := as.integer(curve_id)]
    step_dt[, alpha_ci := conf_int_alpha]

    if (!isTRUE(conf_int)) {
      return(as.data.frame(step_dt))
    }

    by_cols <- if (length(curve_group_cols) > 0) curve_group_cols else NULL

    ci_dt <- data.table::copy(w)
    if (is.null(by_cols)) {
      ci_dt[, xmax := data.table::shift(time, type = "lead")]
    } else {
      ci_dt[, xmax := data.table::shift(time, type = "lead"), by = by_cols]
    }
    ci_dt <- ci_dt[!is.na(xmax)]
    ci_dt[, xmin := time]
    data.table::setnames(ci_dt, "lower", "ymin")
    data.table::setnames(ci_dt, "upper", "ymax")
    ci_dt[, x := NA_real_]
    ci_dt[, y := NA_real_]
    ci_dt[, component := "ci"]
    ci_dt[, group := as.integer(curve_id)]
    ci_dt[, alpha_ci := conf_int_alpha]
    out <- data.table::rbindlist(
      list(step_dt, ci_dt),
      use.names = TRUE,
      fill = TRUE
    )
    print(out)
    as.data.frame(out)
  }
)
