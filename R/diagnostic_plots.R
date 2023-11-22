#' Plot Cook's distances
#'
#' @param model Model with [cooks.distance()] method
#' @param n_highlights Number of points with highest cook distance to highlight in the plot
#' @param cut Boolean indicating to show horizontal line of 4 standard deviations
#' above the Cook's distances mean.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... Further arguments for [graphics::plot()]
#'
#' @return Cook's distances invisible
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_cook(fit)
#' plot_cook(fit, n_highlights = 2)
plot_cook <- function(model, n_highlights = 0, cut = FALSE,
                      xlab = "Index", ylab = "Cook's distance", ...) {
  cook <- stats::cooks.distance(model)
  ylim <- grDevices::extendrange(cook)

  plot(seq_along(cook), cook,
    pch = 20, xlab = xlab, ylab = ylab, ylim = ylim,
    ...
  )
  if (n_highlights > 0) {
    highests <- order(cook, decreasing = TRUE)[seq_len(n_highlights)]
    graphics::text(
      x = highests, y = cook[highests], label = highests,
      pos = 3
    )
  }
  if (cut) {
    graphics::abline(h = base::mean(cook) + 4 * stats::sd(cook))
  }

  return(invisible(cook))
}

#' Plot Residuals against Linear Predictor
#'
#' @param model Model with methods [predict()] or [fitted()]
#' @param residual_fn function to compute residuals, default [rstandard()],
#' other functions could be [rstudent()], [residuals()].
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... Extra arguments to `residual_fn`
#'
#' @return invisible list with plot components
#'
#' @details
#' If the model was fitted using the [glm()] function, it will use the [predict()]
#' method, otherwise, it will use the [fitted()] method.
#'
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_res_vs_linear_predictor(fit)
#' plot_res_vs_linear_predictor(fit, residual_fn = rstudent)
#' plot_res_vs_linear_predictor(fit, residual_fn = residuals)
#'
#'
#' glm_fit <- glm(cyl ~ mpg, family = poisson(), data = mtcars)
#'
#' plot_res_vs_linear_predictor(glm_fit)
#' plot_res_vs_linear_predictor(glm_fit, type = "pearson")
plot_res_vs_linear_predictor <- function(model, residual_fn = stats::rstandard,
                                         xlab = "Linear Predictor",
                                         ylab = "Standardized deviance residuals",
                                         ...) {
  linear_predictor <- switch(class(model)[1],
    glm = stats::predict(model, type = "link"),
    stats::fitted(model)
  )
  y <- residual_fn(model, ...)
  plot(x = linear_predictor, y, xlab = xlab, ylab = ylab)

  return(invisible(list(x = linear_predictor, y = y)))
}
