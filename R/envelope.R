#' Generate Simulated Envelope
#'
#' Generates a normal QQ-plot with simulated envelope residuals.
#'
#' Simulates new responses using [stats::simulate()] and refits the model
#' for each vector of new responses using [get_refit()]. The function then computes
#' residuals for each simulation, sorts them, and constructs envelope bands and
#' a median line based on the quantiles of these residuals.
#'
#' @param model A model that has a [stats::simulate()] method and is compatible with [get_refit()].
#' @param residual_fn A function to calculate model residuals. The default is
#'   [stats::rstudent()] for studentized residuals.
#' @param alpha The significance level for constructing the envelope bounds.
#'   Defaults to 0.05.
#' @param n_sim The number of simulations to perform for envelope construction.
#'   Defaults to 100.
#' @param plot.it Logical. Generate envelope plot.
#' @param ... Extra arguments to [get_refit()]
#'
#' @return An object of class `LD_envelope`, which contains the following components:
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
#' envelope(fit)
#' envelope(fit, residual_fn = function(x) residuals.lm(x, type = "deviance"))
#' envelope(fit, residual_fn = rstandard, n_sim = 200)
#' }
#'
#' @seealso \code{\link{get_refit}}, \code{\link{simulate}}, \code{\link{rstudent}}, \code{\link{plot.LD_envelope}}
#'
#' @export
envelope <- function(model, residual_fn = stats::rstudent,
                     alpha = .05, n_sim = 100,
                     plot.it = TRUE, ...) {
  td <- residual_fn(model)
  n <- length(td)
  e <- matrix(NA, n, n_sim)
  y_star <- stats::simulate(model, n_sim)
  for (i in seq_len(n_sim)) {
    y_ <- y_star[[i]]
    try(
      {
        model_refit <- get_refit(model, y_)
        e[, i] <- sort(residual_fn(model_refit))
      },
      silent = TRUE
    )
  }

  if (anyNA(e)) {
    warning("At least one of the refitted models produced an error.")
  }

  es <- apply(e, 1, stats::quantile, probs = c(alpha / 2, .5, 1 - alpha / 2), na.rm = TRUE)

  tdsort <- sort(td)
  oob <- !(tdsort > es[1, ] & tdsort < es[3, ])

  qq <- stats::qqnorm(td, plot.it = FALSE)
  qq_ord <- order(qq$x)
  qq_sort <- qq$x[qq_ord]

  result <- list(
    expected = qq_sort, observed = qq$y[qq_ord], outside = oob,
    lower = es[1, ], med = es[2, ], upper = es[3, ]
  )

  class(result) <- "LD_envelope"

  if (plot.it) {
    plot(result)
    return(invisible(result))
  }

  result
}

#' Envelope Plot
#'
#' Plot LD_envelope
#'
#' @param x LD_envelope object, usually the result of [envelope()]
#' @param colors Vector of length 2, with color for points outside and inside
#'   the envelope band, respectivelly.
#' @param ylab The label for the y-axis.
#' @param xlab The label for the x-axis.
#' @param ... extra arguments passed to [graphics::plot]
#'
#' @export
plot.LD_envelope <- function(x,
                             colors = c("red", "black"),
                             xlab = "Expected quantiles",
                             ylab = "Observed quantiles", ...) {
  graphics::plot(x$expected, x$observed,
    col = ifelse(x$outside, colors[1], colors[2]),
    type = "p", pch = 20, xlab = xlab, ylab = ylab,
    ...
  )
  graphics::lines(x$expected, x$lower, lty = 1)
  graphics::lines(x$expected, x$upper, lty = 1)
  graphics::lines(x$expected, x$med, lty = 2)
}

#' Generate Simulated Envelope Measures
#'
#' @inheritParams envelope
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
  .Deprecated("envelope")
  envelope(model = model, residual_fn = residual_fn, alpha = alpha, n_sim = n_sim, plot.it = FALSE)
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
  .Deprecated("envelope")
  env_meas <- envelope(model = model, plot.it = FALSE)
  plot(env_meas, colors = colors, xlab = xlab, ylab = ylab, ...)
  return(invisible(env_meas))
}
