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
  model_terms <- stats::terms(stats::formula(model))


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

  return(list(coefs = coefs, variances = variances))
}

compute_p_values <- function(simulation_coefs, simulation_vcov, model) {
    standard_coef <- sapply(seq_along(simulation_coefs), function(i) {
      (simulation_coefs[[i]] - coef(model)) / sqrt(diag(simulation_vcov[[i]]))
    })
    colnames(standard_coef) <- names(simulation_coefs)

  p_values <- 1 - pchisq(standard_coef^2, 1)

  return(p_values)
}

compute_p_values_joint <- function(simulation_coefs, simulation_vcov, model) {
  chisq_stat <- sapply(seq_along(simulation_coefs), function(i) {
    t(simulation_coefs[[i]] - coef(model)) %*%
      solve(simulation_vcov[[i]]) %*%
      (simulation_coefs[[i]] - coef(model))
  })

  p_values <- 1 - pchisq(chisq_stat, length(coef(model)))

  return(p_values)
}

# TODO: Create functions that wraps all functions above
