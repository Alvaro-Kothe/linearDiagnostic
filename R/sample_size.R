simulate_coefficients <- function(model, generator = NULL, n_sim = 100, verbose = FALSE) {
  if (is.null(generator)) {
    y_star <- simulate(model, n_sim)
  } else {
    mu <- fitted(model)
    y_star <- replicate(n_sim,
                        generator(n = length(mu), mu = mu),
                        simplify = FALSE
    )
    names(y_star) <- paste0("sim_", seq_along(y_star))
  }
  model_frame <- model.frame(model)
  model_terms <- terms(formula(model))

  coefs <- list()
  variances <- list()

  for (i in seq_along(y_star)) {
    y_ <- y_star[[i]]

    formula_new_response <- as.formula(
      paste(enquote(y_)[2], "~.", collapse = "")
    )

    refit <- update(model, formula. = formula_new_response, data = model_frame)

    coefs[[i]] <- coef(refit)
    variances[[i]] <- vcov(refit)
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
