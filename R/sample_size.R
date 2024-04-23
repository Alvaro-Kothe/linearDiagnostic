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

compute_p_values <- function(statistic) {
  p_values <- 1 - stats::pchisq(statistic^2, 1)
  p_values
}

compute_p_values_joint <- function(coefs, vcov, generator_coef) {
  used_ginv <- FALSE
  vcov_inv <- tryCatch(solve(vcov),
    error = function(e) {
      used_ginv <<- TRUE
      MASS::ginv(vcov)
    }
  )
  dif_nul <- coefs - generator_coef
  chisq_stat <- t(dif_nul) %*% vcov_inv %*% dif_nul

  p_values <- 1 - stats::pchisq(chisq_stat, length(generator_coef))

  list(p_values, used_ginv)
}

#' Generate P-Values Matrix from Model Coefficient Simulations
#'
#' This function generates p-values by simulating coefficients from
#' a given model and computing Wald test statistics for each simulation.
#'
#' Generate new responses with [stats::simulate()] if generator is `NULL`.
#' If `generator` is provided, it replicates the `generator(model)` by `n_sim`.
#'
#' For each new response calls [get_refit()] to generate a new model with the new
#' response. It gets the fixed effects and the variance and covariance matrix with
#' [get_fixef()] and [get_vcov()].
#'
#' The Univariate Wald test is computed from the Wald statistic, and the p-value
#' is obtained from a chi-squared distribution with one d.f.
#' The joint Wad test computes the inverse from the vcov result and computes the
#' Wald statistic using the `test_coefficients`. The p-values are also obtained
#' from a chi-squared distribution with `length(test_coefficients)` d.f.
#'
#' @param model A model compatible with [get_refit()], [get_fixef()] and [get_vcov()] methods.
#' @param generator An optional function with one argument to generate new response vectors.
#'   If NULL, the model's [simulate()] method is used. Otherwise, the generator
#'   function is used to simulate responses.
#' @param n_sim The number of simulations to perform.
#' @param test_coefficients Numeric vector. A vector with values to be used to compute
#'   the test statistic. It should be the coefficients that was used to compute
#'   the fitted values of the response. If `NULL` defaults to model fixed effects
#' @param ... Additional arguments to be passed to [get_refit()].
#'
#' @return An object of class `LD_pvalues` class, which contains the following components:
#' \describe{
#'   \item{test_coefficients}{Vector of coefficients being tested.}
#'   \item{pvalues_matrix}{Matrix of p-values where each column corresponds to a
#'          simulation and each row corresponds to a coefficient.}
#'   \item{pvalues_joint}{Vector containing the joint p-values obtained from each simulation.}
#'   \item{simulation_fixef}{List of fixed effect coefficient estimates from each simulation.}
#'   \item{simulation_vcov}{List of covariance matrices estimated from each simulation.}
#'   \item{converged}{Logical vector indicating whether model refitting converged for each simulation.}
#'   \item{ginv_used}{Logical vector indicating whether generalized inversed was used to compute the joint p-value.}
#'   \item{responses}{Simulated responses used for refitting the model.}
#' }
#' @seealso [plot.LD_pvalues()] for plotting.
#' @examples
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' generator <- function(object) {
#'   rnorm(nobs(object), mean = fitted.values(object), sd = 2)
#' }
#' get_p_values(model, generator = generator, n_sim = 100)
#'
#' @export
get_p_values <- function(model, n_sim = 1000, generator = NULL,
                         test_coefficients = NULL, ...) {
  out <- list()
  if (is.null(test_coefficients)) {
    test_coefficients <- get_fixef(model)
  }
  out$test_coefficients <- test_coefficients
  out$pvalues_matrix <- matrix(NA, nrow = length(test_coefficients), ncol = n_sim)
  rownames(out$pvalues_matrix) <- names(test_coefficients)
  out$pvalues_joint <- numeric(n_sim)
  out$simulation_fixef <- list()
  out$simulation_vcov <- list()
  out$converged <- rep(NA, n_sim)
  out$ginv_used <- logical(n_sim)
  out$responses <- if (is.null(generator)) {
    stats::simulate(model, n_sim)
  } else {
    replicate(n_sim, generator(model))
  }

  for (i in seq_len(n_sim)) {
    y_star <- out$responses[[i]]
    model_refit <- get_refit(model, y_star, ...)
    fef <- get_fixef(model_refit)
    vc <- get_vcov(model_refit)
    try(out$converged[i] <- get_converged(model_refit), silent = TRUE)
    out$simulation_fixef[[i]] <- fef
    out$simulation_vcov[[i]] <- vc

    statistic <- compute_statistic(coefs = fef, vcov = vc, generator_coef = test_coefficients)
    out$pvalues_matrix[, i] <- compute_p_values(statistic = statistic)

    pj <- compute_p_values_joint(
      coefs = fef,
      vcov = vc,
      generator_coef = test_coefficients
    )
    out$pvalues_joint[i] <- pj[[1]]
    out$ginv_used[i] <- pj[[2]]
  }
  if ((ginv_uses <- sum(out$ginv_used)) > 0) {
    warning(
      "Couldn't inverse vcov from ",
      ginv_uses,
      " simulations and used `MASS::ginv` instead"
    )
  }
  if ((sum_not_conv <- sum(!out$converged, na.rm = TRUE)) > 0) {
    warning("Model didn't converge in ", sum_not_conv, " simulations")
  }

  class(out) <- "LD_pvalues"
  out
}

#' Generate P-Values Matrix from Model Coefficient Simulations
#'
#' This function generates a matrix of p-values by simulating coefficients from
#' a given model and computing test statistics for each simulation.
#'
#' @inheritParams get_p_values
#' @param ... Additional arguments to be passed to `get_p_values`.
#'
#' @return A matrix where each column represents the p-values obtained from a simulation.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' generator <- function(object) {
#'   rnorm(nobs(object), mean = fitted.values(object), sd = 2)
#' }
#' get_p_values_matrix(model, generator = generator, n_sim = 100)
#' }
#' @export
get_p_values_matrix <- function(model, n_sim = 1000,
                                test_coefficients = NULL, ...) {
  .Deprecated("get_p_values")
  get_p_values(model, n_sim = n_sim, test_coefficients = test_coefficients, ...)$pvalues_matrix
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
#' @export
get_p_values_joint <- function(model, n_sim = 1000, test_coefficients = NULL, ...) {
  .Deprecated("get_p_values")
  get_p_values(model, n_sim = n_sim, test_coefficients = test_coefficients, ...)$pvalues_joint
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
#' plot_pvalues_ecdf(fit)
#' }
#' @export
plot_pvalues_ecdf <- function(model,
                              which = seq_along(get_fixef(model)),
                              caption = paste("ECDF of", names(get_fixef(model))[which]),
                              plot_uniform = TRUE,
                              uniform_legend = TRUE,
                              ylab = "Empirical cumulative distribution", xlab = "p-value",
                              args_sim = list(), ...,
                              ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive()) {
  .Deprecated("plot.LD_pvalues")
  margs <- list(model = model)
  full_args <- c(margs, args_sim)
  p_values <- do.call(get_p_values, full_args)
  plot.LD_pvalues(p_values,
    which = which, caption = as.list(caption),
    plot_uniform = plot_uniform, uniform_legend = uniform_legend,
    ylab = ylab, xlab = xlab, ask = ask, ...
  )
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
  .Deprecated("plot.LD_pvalues")
  margs <- list(model = model)
  full_args <- c(margs, args_sim)
  p_values <- do.call(get_p_values, full_args)
  plot.LD_pvalues(p_values,
    which = length(p_values$test_coefficients) + 1,
    plot_uniform = plot_uniform, uniform_legend = uniform_legend,
    ylab = ylab, xlab = xlab, ...
  )
}
