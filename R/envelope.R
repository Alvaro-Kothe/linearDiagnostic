#' Generate Simulated Envelope Measures
#'
#' @param model A model [stats::simulate()] and compatible with [get_refit()].
#' @param residual_fn A function to calculate model residuals. The default is
#'   `stats::rstudent` for studentized residuals.
#' @param alpha The significance level for constructing the envelope bounds.
#'   Defaults to 0.05.
#' @param n_sim The number of simulations to perform for envelope construction.
#'   Defaults to 100.
#' @param ... Extra arguments to [get_refit()]
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{expected}{A vector of expected quantiles from a normal distribution.}
#'   \item{observed}{A vector of observed quantiles from the model residuals.}
#'   \item{outside}{A logical vector indicating whether each observation falls
#'   outside the constructed envelope bounds.}
#'   \item{lower}{The lower bounds of the envelope for each observation.}
#'   \item{med}{The median bounds of the envelope for each observation.}
#'   \item{upper}{The upper bounds of the envelope for each observation.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' envelope_measures(fit, residual_fn = function(x) residuals.lm(x, type = "deviance"))
#'
#' envelope_measures(fit, n_sim = 200)
#' envelope_measures(fit, residual_fn = rstandard, n_sim = 200)
#' }
#'
#' @seealso \code{\link{get_refit}}, \code{\link{simulate}}, \code{\link{rstudent}}
#'
#' @export
envelope_measures <- function(model, residual_fn = stats::rstudent,
                              alpha = .05, n_sim = 100, ...) {
  td <- residual_fn(model)
  n <- length(td)
  e <- matrix(0, n, n_sim)
  y_star <- stats::simulate(model, n_sim)
  for (i in seq_len(n_sim)) {
    y_ <- y_star[[i]]
    model_refit <- get_refit(model, y_)
    e[, i] <- sort(residual_fn(model_refit))
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

plot_envelope_base <- function(expected, observed, lower, med, upper, outside,
                               colors,
                               xlab,
                               ylab) {
  graphics::plot(expected, observed,
    col = ifelse(outside, colors[1], colors[2]),
    type = "p", pch = 20, xlab = xlab, ylab = ylab
  )
  graphics::lines(expected, lower, lty = 1)
  graphics::lines(expected, upper, lty = 1)
  graphics::lines(expected, med, lty = 2)
}

#' Plot Simulated Envelope
#'
#' This function generates a plot to visualize the envelope measures.
#' The plot includes expected quantiles versus observed
#' quantiles, along with lower, median, and upper bounds of the constructed
#' envelope.
#'
#' This function calls [envelope_measures()] and plot it's results.
#'
#' @inheritParams envelope_measures
#' @param colors A vector with two strings, the first element is the color of the
#'   elements outside of the envelope, the second is for the elements inside.
#' @param xlab The label for the x-axis. Defaults to "Expected quantiles".
#' @param ylab The label for the y-axis. Defaults to "Observed quantiles".
#' @param ... Additional arguments for [envelope_measures()]
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{expected}{A vector of expected quantiles from a normal distribution.}
#'   \item{observed}{A vector of observed quantiles from the model residuals.}
#'   \item{outside}{A logical vector indicating whether each observation falls
#'   outside the constructed envelope bounds.}
#'   \item{lower}{The lower bounds of the envelope for each observation.}
#'   \item{med}{The median bounds of the envelope for each observation.}
#'   \item{upper}{The upper bounds of the envelope for each observation.}
#' }
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' plot_envelope(model)
#' }
#'
#' @seealso [graphics::plot()]
#'
#' @export
plot_envelope <- function(model,
                          colors = c("red", "black"),
                          xlab = "Expected quantiles",
                          ylab = "Observed quantiles", ...) {
  env_meas <- envelope_measures(model, ...)

  plot_envelope_base(env_meas$expected, env_meas$observed,
    env_meas$lower, env_meas$med, env_meas$upper, env_meas$outside,
    colors = colors, xlab = xlab, ylab = ylab
  )
  return(invisible(env_meas))
}
