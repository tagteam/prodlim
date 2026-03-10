draw_key_prodlim <- function(data, params, size) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  alpha_val <- data$alpha %||% NA_real_
  if (is.na(alpha_val)) alpha_val <- data$alpha_ci %||% 0.2

  fill_val <- data$fill %||% NA_character_
  col_val  <- data$colour %||% data$color %||% "black"

  if (is.na(fill_val) && !is.na(col_val)) fill_val <- col_val
  if (is.na(fill_val)) fill_val <- "grey70"
  is_stacked <- (length(params$cause) == 1 && params$cause == "stacked")
  if (is_stacked) {
    return(
      grid::rectGrob(
        width = 1, height = 1,
        gp = grid::gpar(
          fill = scales::alpha(fill_val, 1),
          col = NA
        )
      )
    )
  }

  grid::grobTree(
    grid::rectGrob(
      width = 1, height = 1,
      gp = grid::gpar(
        fill = scales::alpha(fill_val, alpha_val),
        col = NA
      )
    ),
    grid::segmentsGrob(
      x0 = 0.1, x1 = 0.9, y0 = 0.5, y1 = 0.5,
      gp = grid::gpar(
        col = col_val,
        lwd = (data$linewidth %||% data$size %||% 0.5) * ggplot2::.pt,
        lty = data$linetype %||% 1
      )
    )
  )
}



GeomProdlim <- ggplot2::ggproto(
  "GeomProdlim",
  ggplot2::Geom,
  required_aes = character(0),

  handle_na = function(data, params) {
    step  <- data[data$component == "step",  , drop = FALSE]
    ci    <- data[data$component == "ci",    , drop = FALSE]
    stack <- data[data$component == "stack", , drop = FALSE]

    if (nrow(step) > 0) {
      step <- step[stats::complete.cases(step[, c("x", "y")]), , drop = FALSE]
    }

    if (nrow(ci) > 0) {
      ci <- ci[
        stats::complete.cases(ci[, c("xmin", "xmax", "ymin", "ymax")]),
        ,
        drop = FALSE
      ]
    }

    if (nrow(stack) > 0) {
      stack <- stack[
        stats::complete.cases(stack[, c("xmin", "xmax", "ymin", "ymax")]),
        ,
        drop = FALSE
      ]
    }

    rbind(step, ci, stack)
  },

  default_aes = ggplot2::aes(
    colour = "black",
    fill = NA,
    linewidth = 0.5,
    linetype = 1,
    alpha = NA
  ),
  draw_key = draw_key_prodlim,

  draw_panel = function(data, panel_params, coord,
                        conf_int = TRUE,
                        conf_int_alpha = 0.2,
                        na.rm = FALSE) {

    `%||%` <- function(a, b) if (!is.null(a)) a else b

    step_data  <- data[data$component == "step",  , drop = FALSE]
    ci_data    <- data[data$component == "ci",    , drop = FALSE]
    stack_data <- data[data$component == "stack", , drop = FALSE]

    grobs <- list()
      
    if (nrow(stack_data) > 0) {
      stack_data$colour <- NA
      stack_data$linewidth <- 0
      if (!("alpha" %in% names(stack_data)) || all(is.na(stack_data$alpha))) {
        stack_data$alpha <- 1
      }
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomRect$draw_panel(
        stack_data,
        panel_params = panel_params,
        coord = coord
      )
      return(do.call(grid::grobTree, grobs))
    }

    if (isTRUE(conf_int) && nrow(ci_data) > 0) {
      if (!("fill" %in% names(ci_data))) {
        ci_data$fill <- ci_data$colour %||% ci_data$color
      } else {
        missing_fill <- is.na(ci_data$fill)
        if (any(missing_fill) && "colour" %in% names(ci_data)) {
          ci_data$fill[missing_fill] <- ci_data$colour[missing_fill]
        }
      }

      ci_data$colour <- NA
      ci_data$linewidth <- 0

      if (!("alpha" %in% names(ci_data)) || all(is.na(ci_data$alpha))) {
        ci_data$alpha <- ci_data$alpha_ci %||% conf_int_alpha
      } else {
        ci_data$alpha[is.na(ci_data$alpha)] <- ci_data$alpha_ci[is.na(ci_data$alpha)]
        ci_data$alpha[is.na(ci_data$alpha)] <- conf_int_alpha
      }

      grobs[[length(grobs) + 1L]] <- ggplot2::GeomRect$draw_panel(
        ci_data,
        panel_params = panel_params,
        coord = coord
      )
    }

    if (nrow(step_data) > 0) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomStep$draw_panel(
        step_data,
        panel_params = panel_params,
        coord = coord,
        direction = "hv"
      )
    }

    do.call(grid::grobTree, grobs)
  }
)
