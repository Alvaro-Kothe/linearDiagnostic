#' Generate Simulated Envelope Measures
#'
#' @param model A model with [stats::update()] method
#' @param residual_fn Function to compute residual using `model`
#' @param alpha Confidence Band level
#' @param n_sim Number of simulations
#'
#' @return A list of vectors with the expected, observed quantiles,
#' boolean, indicating inside or outside of band; and
#' lower band values with `alpha/2`, median, upper band values with `1 - alpha/2`.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' envelope_measures(fit, residual_fn = function(x) residuals.lm(x, type = "deviance"))
#'
#' envelope_measures(fit, n_sim = 200)
#' envelope_measures(fit, residual_fn = rstandard, n_sim = 200)
envelope_measures <- function(model, residual_fn = stats::rstudent,
                              alpha = .05, n_sim = 100) {
  td <- residual_fn(model)
  model_frame <- stats::model.frame(model)
  n <- length(td)
  e <- matrix(0, n, n_sim)
  y_star <- stats::simulate(model, n_sim)
  for (i in seq_len(n_sim)) {
    y_ <- y_star[[i]]

    formula_new_response <- stats::as.formula(
      paste(enquote(y_)[2], "~.", collapse = "")
    )

    refit <- stats::update(model, formula. = formula_new_response, data = model_frame)
    e[, i] <- sort(residual_fn(refit))
  }
  es <- apply(e, 1, function(x) stats::quantile(x, c(alpha / 2, .5, 1 - alpha / 2)))

  tdsort <- sort(td)
  oob <- !(tdsort > es[1, ] & tdsort < es[3, ])

  qq <- stats::qqnorm(td, plot.it = FALSE)
  qq_ord <- order(qq$x)
  qq_sort <- qq$x[qq_ord]

  return(list(
    expected = qq_sort, observed = qq$y[qq_ord], outside = oob,
    lower = es[1, ], med = es[2, ], upper = es[3, ]
  ))
}

plot_envelope_base <- function(expected, observed, lower, med, upper,
                               xlab = "Expected quantiles",
                               ylab = "Observed quantiles") {
  graphics::plot(expected, observed, type = "p", pch = 20, xlab = xlab, ylab = ylab)
  graphics::lines(expected, lower, lty = 1)
  graphics::lines(expected, upper, lty = 1)
  graphics::lines(expected, med, lty = 2)
  return(invisible(list(x = expected, y = observed, lower = lower, upper = upper, med = med)))
}

#' Plot Simulated Envelope
#'
#' @param model A model with [update()] method
#' @param xlab x-axis label

#' @param ylab y-axis label
#' @param ... Additional arguments for [envelope_measures()]
#'
#' @return Creates the simulated envelope plot and returns its components
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_envelope(fit)
plot_envelope <- function(model, xlab = "Expected quantiles", ylab = "Observed quantiles", ...) {
  env_meas <- envelope_measures(model, ...)

  envelope_plot <- plot_envelope_base(env_meas$expected, env_meas$observed,
    env_meas$lower, env_meas$med, env_meas$upper,
    xlab = xlab, ylab = ylab
  )


  return(invisible(envelope_plot))
}
