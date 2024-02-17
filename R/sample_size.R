#' Simulate Model Coefficients and Covariance Matrix
#'
#' This function simulates coefficients and the covariance matrix for a fitted
#' model. It allows for generating new response vectors either using
#' the model's simulation method or a user-specified generator function.
#'
#' @param model A model with [update()], [simulate()], [fitted()],
#' [coef()], [vcov()] methods.
#' @param generator An optional function 2 arguments (n, mu) to generate new response vectors.
#'   If NULL, the model's [simulate()] method is used. Otherwise, the generator
#'   function is used to simulate responses.
#' @param ... Extra arguments to [stats::update()]
#'
#' @return A list containing the simulated coefficients and covariance matrices.
#'   The list includes components:
#' \describe{
#'   \item{coefs}{Coefficient vectors.}
#'   \item{vcov}{Covariance matrices.}
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
simulate_coefficients <- function(model, generator = NULL, ...) {
  if (is.null(generator)) {
    y_star <- stats::simulate(model)[[1]]
  } else {
    mu <- stats::fitted(model)
    y_star <- generator(n = length(mu), mu = mu)
  }

  formula_new_response <- stats::as.formula(
    paste(enquote(y_star)[2], "~.", collapse = "")
  )

  refit <- stats::update(model, formula. = formula_new_response, ...)

  list(coefs = stats::coef(refit), vcov = stats::vcov(refit))
}

#' Compute Wald Statistic
#'
#' Compute Wald statistic for one coefficient using the simulation coefficient
#' and standard error.
#'
#' @param coefs coefficient vector.
#' @param vcov Covariance Matrix, should be a square matrix with p rows.
#' @param generator_coef coefficient vector used to generate `coefs` and `vcov`.
#'
#' @return Wald statistic for each coefficient.
compute_statistic <- function(coefs, vcov, generator_coef) {
  (coefs - generator_coef) / sqrt(diag(vcov))
}

compute_p_values <- function(statistic, df = NULL) {
  if (is.null(df)) {
    p_values <- 1 - stats::pchisq(statistic^2, 1)
  } else {
    p_values <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  }

  p_values
}

compute_p_values_joint <- function(coefs, vcov, generator_coef) {
  vcov_inv <- tryCatch(solve(vcov),
    error = function(e) {
      warning("Couldn't inverse vcov and used `MASS::ginv` instead\n")
      MASS::ginv(vcov)
    }
  )
  dif_nul <- coefs - generator_coef
  chisq_stat <- t(dif_nul) %*% vcov_inv %*% dif_nul

  p_values <- 1 - stats::pchisq(chisq_stat, length(generator_coef))

  return(p_values)
}

#' Generate P-Values Matrix from Model Coefficient Simulations
#'
#' This function generates a matrix of p-values by simulating coefficients from
#' a given model and computing test statistics for each simulation.
#'
#' @inheritParams simulate_coefficients
#' @param n_sim The number of simulations to perform.
#' @param use_tstat Logical. If TRUE, the t-statistic is used for computing p-values.
#'   If FALSE, the z-statistic is used. If NULL, the default is determined based on
#'   the presence of the "t value" column in the summary of the model coefficients.
#' @param ... Additional arguments to be passed to `simulate_coefficients`.
#'
#' @return A matrix where each column represents the p-values obtained from a simulation.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' generator <- function(n, mu) {
#'   rt(n, 3) + mu
#' }
#' get_p_values_matrix(model, generator = generator, n_sim = 100)
#' }
#'
#' @export
get_p_values_matrix <- function(model, n_sim = 1000, use_tstat = NULL, ...) {
  generator_coefs <- stats::coef(model)
  result <- matrix(NA, nrow = length(generator_coefs), ncol = n_sim)
  if (is.null(use_tstat)) {
    use_tstat <- colnames(summary(model)$coefficients)[3] == "t value"
  }

  df <- if (use_tstat) model$df.residual else NULL
  for (i in seq_len(n_sim)) {
    simulation <- simulate_coefficients(
      model = model,
      ...
    )
    statistic <- compute_statistic(
      coefs = simulation$coefs,
      vcov = simulation$vcov,
      generator_coef = generator_coefs
    )
    result[, i] <- compute_p_values(statistic = statistic, df = df)
  }
  result
}

#' Generate P-Values from Model Coefficient Simulations
#'
#' This function generates p-values by simulating coefficients from
#' a given model and computing Wald statistics for each simulation.
#'
#' @inheritParams get_p_values_matrix
#'
#' @return A numeric vector containing the joint p-values obtained from each simulation.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' get_p_values_joint(model, n_sim = 100)
#' }
#'
#' @export
get_p_values_joint <- function(model, n_sim = 1000, ...) {
  generator_coefs <- stats::coef(model)
  result <- double(n_sim)
  ginv_uses <- 0
  for (i in seq_len(n_sim)) {
    simulation <- simulate_coefficients(
      model = model,
      ...
    )
    p_value <- withCallingHandlers(
      compute_p_values_joint(
        coefs = simulation$coefs,
        vcov = simulation$vcov,
        generator_coef = generator_coefs
      ),
      warning = function(warning) {
        ginv_uses <<- ginv_uses + 1
      }
    )
    result[i] <- p_value
  }
  if (ginv_uses > 0) {
    warning(
      "Couldn't inverse vcov from ", ginv_uses, " simulations and used ginv instead"
    )
  }
  result
}

#' Plot Empirical Cumulative Distribution Function (ECDF) of p-values
#'
#' This function generates a series of plots displaying the empirical cumulative
#' distribution function (ECDF) of p-values for selected coefficients in a
#' model. The p-values are computed based on simulated coefficients
#' and covariance matrices.
#'
#' @inheritParams get_p_values_matrix
#' @param which A vector specifying the indices of coefficients to plot. Defaults
#'   to all coefficients.
#' @param caption A character vector providing plot captions for each coefficient.
#' @param plot_uniform Logical. If TRUE, plot uniform distribution.
#' @param uniform_legend Logical. If TRUE, a legend is added to the plot to
#'   distinguish between the p-value and U(0, 1) curves. Defaults to TRUE.
#' @param ylab The label for the y-axis. Defaults to "Empirical cumulative distribution".
#' @param xlab The label for the x-axis. Defaults to "p-value".
#' @param args_sim Extra arguments passed to [get_p_values_matrix()].
#' @param ... Extra arguments from [plot()]
#' @param ask Logical. If TRUE, the user is prompted before each plot. Defaults
#'   to TRUE if in an interactive session and the number of plots is greater
#'   than the available space; otherwise, FALSE.
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
plot_pvalues_ecdf <- function(model,
                              which = seq_along(stats::coef(model)),
                              caption = paste("ECDF of", names(stats::coef(model))[which]),
                              plot_uniform = TRUE,
                              uniform_legend = TRUE,
                              ylab = "Empirical cumulative distribution", xlab = "p-value",
                              args_sim = list(), ...,
                              ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive()) {
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }

  get_p_values_main_args <- list(model = model)
  get_p_values_args <- c(get_p_values_main_args, args_sim)
  p_values <- do.call(get_p_values_matrix, get_p_values_args)
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
#' @param args_sim Extra arguments passed to [get_p_values_joint()].
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
                                    plot_uniform = TRUE,
                                    uniform_legend = TRUE,
                                    ylab = "Empirical cumulative distribution", xlab = "p-value",
                                    args_sim = list(), ...) {
  get_p_values_main_args <- list(model = model)
  get_p_values_args <- c(get_p_values_main_args, args_sim)
  p_value <- do.call(get_p_values_joint, get_p_values_args)
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
