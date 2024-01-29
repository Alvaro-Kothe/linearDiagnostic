#' Simulate Model Coefficients and Covariance Matrix
#'
#' This function simulates coefficients and the covariance matrix for a fitted
#' model. It allows for generating new response vectors either using
#' the model's simulation method or a user-specified generator function. For each
#' simulated response, the model is refitted, and the coefficients and covariance
#' matrix of the refitted model are stored.
#'
#' @param model A model with [update()], [simulate()], [fitted()],
#' [coef()], [vcov()] methods.
#' @param generator An optional function 2 arguments (n, mu) to generate new response vectors.
#'   If NULL, the model's [simulate()] method is used. Otherwise, the generator
#'   function is used to simulate responses.
#' @param n_sim The number of simulations to perform. Defaults to 100.
#' @param ... Extra arguments to [stats::update()]
#'
#' @return A list containing the simulated coefficients and covariance matrices.
#'   The list includes components:
#' \describe{
#'   \item{coefs}{A list of coefficient vectors for each simulated response.}
#'   \item{vcov}{A list of covariance matrices for each simulated response.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lm(mpg ~ cyl, data = mtcars)
#' generator <- function(n, mu) rt(n, 3) + mu
#' simulate_coefficients(fit, generator = generator)
#' }
#'
#' @export
simulate_coefficients <- function(model, generator = NULL, n_sim = 100, ...) {
  if (is.null(generator)) {
    y_star <- stats::simulate(model, n_sim)
  } else {
    mu <- stats::fitted(model)
    y_star <- replicate(n_sim,
      generator(n = length(mu), mu = mu),
      simplify = FALSE
    )
    names(y_star) <- paste0("sim_", seq_along(y_star))
  }

  coefs <- list()
  variances <- list()

  for (i in seq_along(y_star)) {
    y_ <- y_star[[i]]

    formula_new_response <- stats::as.formula(
      paste(enquote(y_)[2], "~.", collapse = "")
    )

    refit <- stats::update(model, formula. = formula_new_response, ...)

    coefs[[i]] <- stats::coef(refit)
    variances[[i]] <- stats::vcov(refit)
  }

  names(coefs) <- names(variances) <- names(y_star)

  return(list(coefs = coefs, vcov = variances))
}

compute_statistic <- function(simulation_coefs, simulation_vcov, generator_coef) {
  standard_coef <- sapply(seq_along(simulation_coefs), function(i) {
    (simulation_coefs[[i]] - generator_coef) / sqrt(diag(simulation_vcov[[i]]))
  })
  colnames(standard_coef) <- names(simulation_coefs)

  return(standard_coef)
}

compute_p_values <- function(statistic, df = NULL) {
  if (is.null(df)) {
    p_values <- 1 - stats::pchisq(statistic^2, 1)
  } else {
    p_values <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  }

  return(p_values)
}

compute_p_values_joint <- function(simulation_coefs, simulation_vcov, generator_coef) {
  ginv_uses <- 0
  chisq_stat <- sapply(seq_along(simulation_coefs), function(i) {
    vcov_inv <- tryCatch(solve(simulation_vcov[[i]]),
      error = function(e) {
        ginv_uses <<- ginv_uses + 1
        MASS::ginv(simulation_vcov[[i]])
      }
    )
    dif_nul <- simulation_coefs[[i]] - generator_coef
    t(dif_nul) %*%
      vcov_inv %*%
      dif_nul
  })

  if (ginv_uses > 0) {
    warning(
      "Couldn't inverse vcov from ",
      ginv_uses,
      " simulations and used ginv instead\n"
    )
  }
  p_values <- 1 - stats::pchisq(chisq_stat, length(generator_coef))

  return(p_values)
}

#' Plot Empirical Cumulative Distribution Function (ECDF) of p-values
#'
#' This function generates a series of plots displaying the empirical cumulative
#' distribution function (ECDF) of p-values for selected coefficients in a
#' model. The p-values are computed based on simulated coefficients
#' and covariance matrices.
#'
#' @inheritParams simulate_coefficients
#' @param which A vector specifying the indices of coefficients to plot. Defaults
#'   to all coefficients.
#' @param caption A character vector providing plot captions for each coefficient.
#' @param plot_uniform Logical. If TRUE, plot uniform distribution.
#' @param uniform_legend Logical. If TRUE, a legend is added to the plot to
#'   distinguish between the p-value and U(0, 1) curves. Defaults to TRUE.
#' @param ylab The label for the y-axis. Defaults to "Empirical cumulative distribution".
#' @param xlab The label for the x-axis. Defaults to "p-value".
#' @param args_sim Extra arguments passed to [simulate_coefficients()].
#' @param ... Extra arguments from [plot()]
#' @param ask Logical. If TRUE, the user is prompted before each plot. Defaults
#'   to TRUE if in an interactive session and the number of plots is greater
#'   than the available space; otherwise, FALSE.
#' @param use_tstat Logical. If TRUE, the t-statistic is used for computing p-values.
#'   If FALSE, the z-statistic is used. If NULL, the default is determined based on
#'   the presence of the "t value" column in the summary of the model coefficients.
#'
#' @return Matrix with p_values obtained from the simulations.
#'
#' @examples
#' \dontrun{
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_pvalues_ecdf(fit)
#' }
#' @export
plot_pvalues_ecdf <- function(model, generator = NULL, n_sim = 1000,
                              which = seq_along(stats::coef(model)),
                              caption = paste("ECDF of", names(stats::coef(model))[which]),
                              plot_uniform = TRUE,
                              uniform_legend = TRUE,
                              ylab = "Empirical cumulative distribution", xlab = "p-value",
                              args_sim = list(), ...,
                              ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive(),
                              use_tstat = NULL) {
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  simulation <- do.call(simulate_coefficients, c(
    model = list(model), generator = list(generator),
    n_sim = n_sim, args_sim
  ))
  statistic <- compute_statistic(
    simulation_coefs = simulation$coefs,
    simulation_vcov = simulation$vcov,
    generator_coef = stats::coef(model)
  )

  if (is.null(use_tstat)) {
    use_tstat <- colnames(summary(model)$coefficients)[3] == "t value"
  }
  df_ <- if (use_tstat) model$df.residual else NULL
  p_values <- compute_p_values(statistic = statistic, df = df_)
  for (coef_idx in seq_along(which)) {
    i <- which[coef_idx]
    p_value <- p_values[i, ]
    ecdf_ <- stats::ecdf(p_value)
    x <- seq(-0.01, 1.01, length.out = 201)
    plot(x, ecdf_(x), type = "l", main = caption[coef_idx], ylab = ylab, xlab = xlab, ...)
    if (plot_uniform) {
      graphics::lines(x, stats::punif(x), lty = 2, col = "gray30")
      if (uniform_legend) {
        graphics::legend("topleft",
          legend = c("p-value", "U(0, 1)"),
          lty = c(1, 2),
          col = c("black", "gray30")
        )
      }
    }
  }
  return(invisible(p_values))
}

#' Plot Joint Empirical Cumulative Distribution Function (ECDF) of p-values
#'
#' This function generates a plot displaying the joint empirical cumulative
#' distribution function (ECDF) of p-values for all coefficients in a model.
#' The p-values are computed based on simulated coefficients and covariance
#' matrices.
#'
#' @inheritParams plot_pvalues_ecdf
#'
#' @return A vector of joint p-values for all coefficients.
#'
#' @examples
#' \dontrun{
#' fit <- lm(mpg ~ cyl, data = mtcars)
#' plot_joint_pvalues_ecdf(fit)
#' }
#' @export
plot_joint_pvalues_ecdf <- function(model,
                                    generator = NULL, n_sim = 1000,
                                    plot_uniform = TRUE,
                                    uniform_legend = TRUE,
                                    ylab = "Empirical cumulative distribution", xlab = "p-value",
                                    args_sim = list(), ...) {
  simulation <- do.call(simulate_coefficients, c(
    model = list(model), generator = list(generator),
    n_sim = n_sim, args_sim
  ))
  p_value <- compute_p_values_joint(
    simulation_coefs = simulation$coefs,
    simulation_vcov = simulation$vcov,
    generator_coef = stats::coef(model)
  )
  ecdf_ <- stats::ecdf(p_value)
  x <- seq(-0.01, 1.01, length.out = 201)
  plot(x, ecdf_(x), type = "l", ylab = ylab, xlab = xlab, ...)
  if (plot_uniform) {
    graphics::lines(x, stats::punif(x), lty = 2, col = "gray30")
    if (uniform_legend) {
      graphics::legend("topleft",
        legend = c("p-value", "U(0, 1)"),
        lty = c(1, 2),
        col = c("black", "gray30")
      )
    }
  }
  return(invisible(p_value))
}
