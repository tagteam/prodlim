stat_prodlim_setup_data <- function(data, params) {
    data <- data.table::as.data.table(data)
    data
}

stat_prodlim_compute_panel <- function(data, scales,
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
      "linewidth",
      "component", "Event", "stacked", "alpha_ci"
    )

    model_cols <- setdiff(names(data), ignore_formula)

    ff <- stats::formula(
      paste0(
        "prodlim::Hist(x, event, cens.code = '", cens.code, "') ~ ",
        if (length(model_cols) == 0) "1" else paste(model_cols, collapse = "+")
      )
    )

    fit <- prodlim::prodlim(formula = ff, data = data, type = type)

    summary_args <- list(
      x = fit,
      surv = identical(type, "surv"),
      percent = percent
    )

    stacked <- identical(cause, "stacked")
    multi_cause <- FALSE
    if (attr(fit$model.response, "model") == "competing.risks") {
        if (stacked || length(cause) == 0) {
            if (stacked){
                requested_causes <- attr(fit$model.response, "states")
            } else{
                requested_causes <- attr(fit$model.response, "states")[[1]]
            }
        } else {
            cause <- intersect(cause, attr(fit$model.response, "states"))
            if (length(cause) == 0) {
                stop(paste0("Cannot find cause(s) ", paste0(cause,collapse = ", "),". The values of data$event are ",paste0(attr(fit$model.response, "states"), collapse = ", ")))
            }
            requested_causes <- cause
        }
        summary_args$cause <- requested_causes
    } else {
        requested_causes <- 1
        stacked <- FALSE
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
            data.table::setorderv(w, c(if (length(model_cols)) model_cols else NULL, "cause", "time"))
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

    by_cols <- if (isTRUE(multi_cause)) c(model_cols, "cause") else model_cols
    data.table::setorderv(w, c(by_cols, "time"))

    if (length(by_cols) == 0) {
        data.table::set(w, j = "Event", value = factor(rep("1", nrow(w))))
    } else {
        if (identical(by_cols[[1]], "cause")) {
            data.table::set(w, j = "Event", value = factor(w[["cause"]], levels = requested_causes))
        } else {
            curve_labels <- as.character(w[[by_cols[[1]]]])
            if (length(by_cols) > 1) {
                for (j in 2:length(by_cols)) {
                    curve_labels <- paste(curve_labels, w[[by_cols[[j]]]], sep = ", ")
                }
            }
            data.table::set(w, j = "Event", value = factor(curve_labels))
        }
    }

    if (isTRUE(stacked) && isTRUE(multi_cause)) {
        if (length(model_cols) == 0) {
            stack_dt <- data.table::copy(w)
            data.table::setorder(stack_dt, cause, time)
            stack_dt[, xmax := data.table::shift(time, type = "lead"),by = cause]
            data.table::setorder(stack_dt, time, cause)
            stack_dt[, ymin := data.table::shift(cumsum(absolute_risk), fill=0), by=time]
            stack_dt[, ymax := cumsum(absolute_risk), by=time]
            data.table::setnames(stack_dt, "time", "xmin")
        } else {
            stop("Cannot split stacked plot.")
        }
        stack_dt <- stack_dt[!is.na(xmax)]
        stack_dt[, x := NA_real_]
        stack_dt[, y := NA_real_]
        stack_dt[, component := "stack"]
        stack_dt[, group := as.integer(Event)]
        stack_dt[, alpha_ci := 1]
        return(as.data.frame(stack_dt))
    }
    step_dt <- data.table::copy(w)
    data.table::setnames(step_dt, y_name, "y")
    data.table::setnames(step_dt, "time", "x")
    step_dt[, component := "step"]
    step_dt[, xmin := NA_real_]
    step_dt[, xmax := NA_real_]
    step_dt[, ymin := NA_real_]
    step_dt[, ymax := NA_real_]
    step_dt[, group := as.integer(Event)]
    step_dt[, alpha_ci := conf_int_alpha]

    if (!isTRUE(conf_int)) {
      return(as.data.frame(step_dt))
    }

    ci_dt <- data.table::copy(w)
    if (length(model_cols) == 0) {
      ci_dt[, xmax := data.table::shift(time, type = "lead"), by = "cause"]
    } else {
      ci_dt[, xmax := data.table::shift(time, type = "lead"), by = c(model_cols, "cause")]
    }
    ci_dt <- ci_dt[!is.na(xmax)]
    ci_dt[, xmin := time]
    data.table::setnames(ci_dt, "lower", "ymin")
    data.table::setnames(ci_dt, "upper", "ymax")
    ci_dt[, x := NA_real_]
    ci_dt[, y := NA_real_]
    ci_dt[, component := "ci"]
    ci_dt[, group := as.integer(Event)]
    ci_dt[, alpha_ci := conf_int_alpha]

    out <- data.table::rbindlist(
      list(step_dt, ci_dt),
      use.names = TRUE,
      fill = TRUE
    )
    as.data.frame(out)
}

StatProdlim <- ggplot2::ggproto(
    "StatProdlim",
    ggplot2::Stat,
    required_aes = c("x", "event"),
    extra_params = c("na.rm",
                     "type",
                     "cause",
                     "conf_int",
                     "percent",
                     "times",
                     "timeconverter",
                     "cens.code",
                     "conf_int_alpha"
    ),
    setup_data = stat_prodlim_setup_data,
    compute_panel = stat_prodlim_compute_panel
)
