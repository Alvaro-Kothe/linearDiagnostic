#' Simulate Model Coefficients
#'
#' @param model A model with [update()] method
#' @param generator Optional custom generator function with 2 arguments (n, mu)
#' @param n_sim Number of simulations
#'
#' @return A list with a list for coefficients and a list with vcov, each of length n_sim
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' generator <- function(n, mu) rt(n, 3) + mu
#'
#' simulate_coefficients(fit, generator = generator)
#'
simulate_coefficients <- function(model, generator = NULL, n_sim = 100) {
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
  model_frame <- stats::model.frame(model)

  coefs <- list()
  variances <- list()

  for (i in seq_along(y_star)) {
    y_ <- y_star[[i]]

    formula_new_response <- stats::as.formula(
      paste(enquote(y_)[2], "~.", collapse = "")
    )

    refit <- stats::update(model, formula. = formula_new_response, data = model_frame)

    coefs[[i]] <- stats::coef(refit)
    variances[[i]] <- stats::vcov(refit)
  }

  names(coefs) <- names(variances) <- names(y_star)

  return(list(coefs = coefs, vcov = variances))
}

compute_p_values <- function(simulation_coefs, simulation_vcov, generator_coef) {
  standard_coef <- sapply(seq_along(simulation_coefs), function(i) {
    (simulation_coefs[[i]] - generator_coef) / sqrt(diag(simulation_vcov[[i]]))
  })
  colnames(standard_coef) <- names(simulation_coefs)

  p_values <- 1 - stats::pchisq(standard_coef^2, 1)

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

#' Plot Simulation p-values
#'
#' @param model A model with [update()] method.
#' @param generator Optional custom generator function with 2 arguments (n, mu).
#' @param n_sim Number of simulations.
#' @param which vector of integers with model coefficients indices to plot.
#' @param caption Plot title for each coefficient.
#' @param ylab y-axis label
#' @param xlab x-axis label
#' @param ... Extra arguments from [plot()]
#' @param ask Logical, ask to show next plot
#'
#' @return Matrix with p_values obtained from the simulations.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_pvalues_ecdf(fit)
plot_pvalues_ecdf <- function(model, generator = NULL, n_sim = 1000,
                              which = seq_along(stats::coef(model)),
                              caption = paste("ECDF of", names(stats::coef(model)))[which],
                              ylab = "Fn(x)", xlab = "x", ...,
                              ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive()) {
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  simulation <- simulate_coefficients(model = model, generator = generator, n_sim = n_sim)
  p_values <- compute_p_values(
    simulation_coefs = simulation$coefs,
    simulation_vcov = simulation$vcov,
    generator_coef = stats::coef(model)
  )
  for (i in which) {
    p_value <- p_values[i, ]
    ecdf_ <- stats::ecdf(p_value)
    x <- seq(-0.01, 1.01, length.out = 201)
    plot(x, ecdf_(x), type = "l", main = caption[i], ylab = ylab, xlab = xlab, ...)
    graphics::lines(x, stats::punif(x), lty = 2)
  }
  return(invisible(p_values))
}

#' Plot Simulation joint p-values
#'
#' @param model A model with [update()] method.
#' @param generator Optional custom generator function with 2 arguments (n, mu).
#' @param n_sim Number of simulations.
#' @param ylab y-axis label
#' @param xlab x-axis label
#' @param ... Extra arguments from [plot()]
#'
#' @return Double vector with joint p-values.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#'
#' plot_joint_pvalues_ecdf(fit)
plot_joint_pvalues_ecdf <- function(model, generator = NULL, n_sim = 1000,
                                    ylab = "Fn(x)", xlab = "x", ...) {
  simulation <- simulate_coefficients(model = model, generator = generator, n_sim = n_sim)
  p_value <- compute_p_values_joint(
    simulation_coefs = simulation$coefs,
    simulation_vcov = simulation$vcov,
    generator_coef = stats::coef(model)
  )
  ecdf_ <- stats::ecdf(p_value)
  x <- seq(-0.01, 1.01, length.out = 201)
  plot(x, ecdf_(x), type = "l", ylab = ylab, xlab = xlab, ...)
  graphics::lines(x, stats::punif(x), lty = 2)
  return(invisible(p_value))
}
