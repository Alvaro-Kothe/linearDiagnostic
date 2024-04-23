#' Plot Cook's distances
#'
#' @param model Model with [cooks.distance()] method
#' @param n_highlights The number of observations with the highest Cook's
#'   distance to highlight on the plot. Defaults to 0 (no highlights).
#' @param cut Logical. If TRUE, adds a cutoff line at the mean plus four times
#'   the standard deviation of Cook's distance. Defaults to FALSE.
#' @param xlab The label for the x-axis. Defaults to "Index".
#' @param ylab The label for the y-axis. Defaults to "Cook's distance".
#' @param ... Further arguments for [graphics::plot()]
#'
#' @return An invisible object representing Cook's distance values.
#'
#' @examples
#' \dontrun{
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_cook(fit)
#' plot_cook(fit, n_highlights = 2)
#' }
#'
#' @seealso \code{\link{cooks.distance}}, \code{\link{plot}}
#' @export
plot_cook <- function(model, n_highlights = 0, cut = FALSE,
                      xlab = "Index", ylab = "Cook's distance", ...) {
  cook <- stats::cooks.distance(model)
  ylim <- grDevices::extendrange(cook)

  plot(seq_along(cook), cook,
    pch = 20, xlab = xlab, ylab = ylab, ylim = ylim,
    ...
  )
  if (n_highlights > 0) {
    highests <- order(cook, decreasing = TRUE)[seq_len(min(c(n_highlights, length(cook))))]
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
#' @param residual_fn A function to calculate model residuals. The default is
#'   `stats::rstandard`.
#' @param xlab The label for the x-axis. Defaults to "Linear Predictor".
#' @param ylab The label for the y-axis. Defaults to "Standardized deviance residuals".
#' @param ... Extra arguments to `residual_fn` and [plot()].
#'
#' @return An invisible list containing the linear predictor (x) and standardized
#'   deviance residuals (y).
#'
#' @details
#' If the model was fitted using the [glm()] function, it will use the [predict()]
#' method with `type = link`, otherwise, it will use the [fitted()] method.
#'
#' @examples
#' \dontrun{
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
#' }
#'
#' @export
plot_res_vs_linear_predictor <- function(model, residual_fn = stats::rstandard,
                                         xlab = "Linear Predictor",
                                         ylab = "Standardized deviance residuals",
                                         ...) {
  linear_predictor <- switch(class(model)[1],
    glm = stats::predict(model, type = "link"),
    stats::fitted(model)
  )
  y <- residual_fn(model, ...)
  plot(x = linear_predictor, y, xlab = xlab, ylab = ylab, ...)

  return(invisible(list(x = linear_predictor, y = y)))
}

#' Plot Empirical Cumulative Distribution Function (ECDF) of p-values
#'
#' This function creates several plots with the empirical cumulative distribution
#' of the p-values obtained thorugh simulation.
#'
#' If the asymptotic approximation is valid the distribution of the p-values
#' should be close to an uniform distribution.
#' Discrepancies are highlighted, by default it verifies the significance on the
#' most commonly used significance values are 0.01, 0.05 and 0.10.
#'
#' The reported KS (Kolmogorov-Smirnov) test is the result of the "two-sided" [stats::ks.test()] function
#' comparing the observed p-values distribution with the uniform.
#' The test may reject the KS test due to few simulations, make sure that the lines
#' shown in the plot are smooth before drawing any conclusions.
#'
#' @param x LD_pvalues object, usually the result of [get_p_values()]
#' @param which A vector specifying the indices of coefficients to plot.
#'  If index is bigger than the number of coefficients it plots the joint p_value.
#' @param caption A list with caption for each plot.
#' @param ks_test If `TRUE` inserts Kolmogorov-Smirnov p-value in the graphic.
#' @param signif Points to verify discrepancy.
#' @param discrepancy_tol Threshold to consider point discrepant.
#' @param plot_uniform Logical. If TRUE, plot uniform distribution.
#' @param uniform_legend Logical. If TRUE, a legend is added to the plot to
#'   distinguish between the p-value and U(0, 1) curves. Defaults to TRUE.
#' @param ylab The label for the y-axis. Defaults to "Empirical cumulative distribution".
#' @param xlab The label for the x-axis. Defaults to "p-value".
#' @param ... extra arguments passed to [graphics::plot]
#' @param ask Logical. If TRUE, the user is prompted before each plot. Defaults
#'   to TRUE if in an interactive session and the number of plots is greater
#'   than the available space; otherwise, FALSE.
#'
#' @return A vector of joint p-values for all coefficients.
#'
#' @examples
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' p_values_ld <- get_p_values(model, n_sim = 100)
#' plot(p_values_ld)
#' @export
plot.LD_pvalues <- function(x,
                            which = seq_len(length(x$test_coefficients) + 1),
                            caption = as.list(paste("ECDF of", c(names(x$test_coefficients), "all coefficients"))),
                            ks_test = TRUE, signif = c(0.01, 0.05, 0.10),
                            discrepancy_tol = .1,
                            plot_uniform = TRUE, uniform_legend = TRUE,
                            ylab = "Empirical cumulative distribution", xlab = "p-value",
                            ...,
                            ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive()) {
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  for (i in which) {
    p_values <- if (i > length(x$test_coefficients)) x$pvalues_joint else x$pvalues_matrix[i, ]
    plot_ecdf_pvalue(p_values,
      ks_test = ks_test, signif = signif,
      discrepancy_tol = discrepancy_tol,
      plot_uniform = plot_uniform, uniform_legend = uniform_legend,
      main = caption[[i]], xlab = xlab, ylab = ylab, ...
    )
  }
}

#' Plot Empirical Cumulative Distribution Function (ECDF) of p-values
#'
#' @inheritParams plot.LD_pvalues
#' @param p_values vector of p-values
#' @param main main caption passed to [plot]
#' @export
plot_ecdf_pvalue <- function(p_values,
                             ks_test, signif,
                             discrepancy_tol,
                             plot_uniform, uniform_legend,
                             main, ylab, xlab, ...) {
  ecdf_ <- stats::ecdf(p_values)
  alpha_ <- seq(-0.01, 1.01, length.out = 201)
  plot(
    alpha_, ecdf_(alpha_),
    type = "l",
    main = main, ylab = ylab, xlab = xlab, ...
  )
  if (ks_test) {
    ks_res <- test_uniform_dist(p_values)$p.value
    graphics::legend("top", legend = sprintf("KS p-value: %.5f", ks_res), bty = "n")
  }
  discrepancies <- character()
  for (expected_rejection in signif) {
    observed <- ecdf_(expected_rejection)
    discrepancy <- (observed / expected_rejection) - 1.0
    if (abs(discrepancy) > discrepancy_tol) {
      discrepancies <- c(discrepancies, sprintf(
        "%.2f: %.3f (%+.0f%%)",
        expected_rejection, observed, discrepancy * 100
      ))
      graphics::segments(expected_rejection, 0.0, y1 = observed, lty = 2)
    }
  }
  if (length(discrepancies) > 0) {
    graphics::legend("bottomright", discrepancies, bty = "n")
  }
  if (plot_uniform) {
    graphics::lines(alpha_, stats::punif(alpha_), lty = 2, col = "gray30")
    if (uniform_legend) {
      graphics::legend("topleft",
        legend = c("p-value", "U(0, 1)"),
        lty = c(1, 2),
        col = c("black", "gray30")
      )
    }
  }
}


test_uniform_dist <- function(x) {
  stats::ks.test(x, stats::punif)
}
