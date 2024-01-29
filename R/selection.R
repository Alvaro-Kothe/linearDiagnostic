remove_variable <- function(x) {
  intercept <- "(Intercept)"
  remove_intercept <- intercept %in% x

  rhs <- paste0(x, collapse = "-")
  if (remove_intercept) rhs <- paste(rhs, "+ 0")
  out <- stats::as.formula(paste0(". ~ . - ", rhs))

  return(out)
}

add_variable <- function(x) {
  rhs <- paste0(x, collapse = "+")
  out <- stats::as.formula(paste0(". ~ . + ", rhs))

  return(out)
}

#' Select covariates backwise
#'
#' @inheritParams select_covariates
#'
#' @return the selected model if `return_step_results` is `False`.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ ., data = mtcars)
#' backward_selection(model)
#' }
#'
#' @export
backward_selection <- function(model, threshold = .15,
                               measure_fn = function(x) summary(x)[["coefficients"]][, 4],
                               data = NULL,
                               max_steps = 1000,
                               return_step_results = FALSE,
                               do_not_remove = c("(Intercept)")) {
  return(select_covariates(
    model = model, direction = "backward",
    threshold = threshold, measure_fn = measure_fn, data = data, max_steps = max_steps,
    return_step_results = return_step_results,
    do_not_remove = do_not_remove
  ))
}


#' Select covariates bidirectional
#'
#' @inheritParams select_covariates
#'
#' @return the selected model if `return_step_results` is `False`.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ ., data = mtcars)
#' select_covariates(model)
#'
#' ## measure_fn with two parameters
#' fit <- lm(mpg ~ ., data = mtcars)
#'
#' lrt <- function(model1, model2) {
#'   lrt_stat <- 2 * (logLik(model1)[1L] - logLik(model2)[1L])
#'   return(1 - pchisq(lrt_stat, 1))
#' }
#'
#' select_covariates(fit, measure_fn = lrt)
#' }
#'
#' @export
bidirectional_selection <- function(model, threshold = .15,
                                    addable_coefs = NULL,
                                    measure_fn = function(x) summary(x)[["coefficients"]][, 4],
                                    data = NULL,
                                    max_steps = 1000,
                                    return_step_results = FALSE,
                                    do_not_remove = c("(Intercept)")) {
  return(select_covariates(
    model = model, direction = "both",
    threshold = threshold, addable_coefs = addable_coefs,
    measure_fn = measure_fn, data = data, max_steps = max_steps,
    return_step_results = return_step_results,
    do_not_remove = do_not_remove
  ))
}

backward_values <- function(model,
                            measure_fn,
                            measure_one_at_time,
                            data,
                            do_not_remove,
                            seen_states) {
  nargs_measure_fn <- length(formals(measure_fn))
  cur_coefs <- names(stats::coef(model))
  candidates_remove <- stats::setNames(nm = setdiff(cur_coefs, do_not_remove))

  evaluate_one_at_time <- measure_one_at_time || nargs_measure_fn == 2L

  possible_next_states <- lapply(candidates_remove, function(x) {
    names(stats::coef(model))[names(stats::coef(model)) != x]
  })

  next_state_seen <- sapply(
    possible_next_states,
    function(x) any(sapply(seen_states, setequal, x))
  )

  next_states_lookup <- names(next_state_seen)[!next_state_seen]

  if (!evaluate_one_at_time) {
    values <- measure_fn(model)[next_states_lookup]
  } else {
    values <- sapply(next_states_lookup, function(remove_candidate) {
      next_formula <- remove_variable(remove_candidate)

      next_fit <- stats::update(model, formula. = next_formula, data = data)

      eval_value <- if (nargs_measure_fn == 1L) {
        measure_fn(next_fit)
      } else if (nargs_measure_fn == 2L) {
        measure_fn(model, next_fit)
      }

      return(eval_value)
    })
  }

  return(values)
}

forward_values <- function(model,
                           measure_fn,
                           measure_one_at_time,
                           data,
                           addable_coefs,
                           seen_states) {
  nargs_measure_fn <- length(formals(measure_fn))
  add_candidates <- setdiff(addable_coefs, names(stats::coef(model)))
  possible_next_states <- lapply(add_candidates, c, names(stats::coef(model)))
  next_state_seen <- sapply(
    possible_next_states,
    function(x) any(sapply(seen_states, setequal, x))
  )
  add_candidates <- add_candidates[!next_state_seen]
  evaluate_one_at_time <- measure_one_at_time || nargs_measure_fn == 2L

  values <- sapply(add_candidates, function(add_cand) {
    next_formula <- add_variable(add_cand)

    next_fit <- stats::update(model, formula. = next_formula, data = data)

    if (!evaluate_one_at_time) {
      return(measure_fn(next_fit)[[add_cand]])
    } else {
      if (nargs_measure_fn == 2L) {
        return(measure_fn(next_fit, model))
      } else {
        return(measure_fn(next_fit))
      }
    }
  })

  return(values)
}

update_model_remove <- function(model,
                                values,
                                threshold,
                                data,
                                minimize_only = FALSE) {
  to_remove <- if (minimize_only) which.min(values) else which.max(values)
  if (length(to_remove) == 0L) {
    return(NULL)
  }

  value <- values[to_remove]
  cant_remove <- if (!minimize_only) value <= threshold else value >= threshold
  if (cant_remove) {
    return(NULL)
  }

  removed_var <- names(to_remove)
  next_formula <- remove_variable(removed_var)

  model <- stats::update(model, formula. = next_formula, data = data)

  return(list(fit = model, removed_var = value))
}

update_model_add <- function(model, values, threshold, data) {
  to_add <- which.min(values)
  if (length(to_add) == 0L) {
    return(NULL)
  }

  value <- values[to_add]

  if (value > threshold) {
    return(NULL)
  }

  added_var <- names(to_add)
  next_formula <- add_variable(added_var)

  model <- stats::update(model, formula. = next_formula, data = data)

  return(list(fit = model, added_var = value))
}

#' Select covariates
#'
#' @param model A model with [stats::update()], [stats::coef()] methods.
#' @param threshold Value threshold to remove variable. It can be a fixed value
#' or a function. The variable is removed if `measure_fn(model) > threshold` and
#' added if `measure_fn(model) <= threshold`.
#' @param direction The direction of variable selection. Options include "backward",
#'   "forward", or "both". Defaults to "both".
#' @param addable_coefs A vector of coefficients that can be added during forward selection.
#'   Defaults to all coefficients in the model.
#' @param measure_fn Function with model as argument and returns values to be used by
#' `threshold`. It can also compare two models, where during forward step
#' it calls `measure_fn(candidate_model, current_selected_model)` and
#' during backward step it calls `measure_fn(current_selected_model, candidate_model)`.
#' Defaults to the p-value from the summary of the coefficients.
#' @param measure_one_at_time Boolean indicating to apply `measure_fn` to each
#' variable individually during forward and backward steps.
#' Set this option to `TRUE` if `measure_fn` returns an atomic value, for example if
#' `measure_fn` is `AIC`.
#' @param minimize_only Logical indicating that during backward model update
#' it should minimize the `measure_fn` instead of maximize it.
#' @param data Data to be used for model refit.
#' @param max_steps The maximum number of steps for the variable selection process.
#'   Defaults to 1000.
#' @param return_step_results Logical. If TRUE, the function returns a list
#'   containing the final fitted model and a log of the selection steps.
#'   Defaults to FALSE.
#' @param do_not_remove A character vector specifying variables that should not
#'   be removed during backward selection. Defaults to "(Intercept)".
#'
#' @return A fitted model with selected covariates based on the variable selection process.
#'   If \code{return_step_results} is TRUE, a list containing the final fitted model
#'   and a log of the selection steps is returned.
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ ., data = mtcars)
#' select_covariates(model)
#'
#' ## measure_fn with two parameters
#'
#' lrt <- function(model1, model2) {
#'   lrt_stat <- 2 * (logLik(model1)[1L] - logLik(model2)[1L])
#'   return(1 - pchisq(lrt_stat, 1))
#' }
#'
#' select_covariates(model, measure_fn = lrt)
#'
#' ## AICc selection
#'
#' AICc <- function(model) {
#'   loglike <- logLik(model)
#'   df <- attr(loglike, "df")
#'   nobs <- attr(loglike, "nobs")
#'   aic <- -2 * as.numeric(loglike) + 2 * df
#'
#'   aicc <- aic + (2 * (df^2) + 2 * df) / (nobs - df - 1)
#'
#'   return(aicc)
#' }
#'
#' selection <- select_covariates(model,
#'   measure_fn = AICc,
#'   threshold = AICc,
#'   measure_one_at_time = TRUE,
#'   minimize_only = TRUE,
#'   direction = "both",
#'   data = mtcars
#' )
#' }
#'
#' @export
select_covariates <- function(model,
                              threshold = .15,
                              direction = "both",
                              addable_coefs = names(stats::coef(model)),
                              measure_fn = function(x) summary(x)[["coefficients"]][, 4],
                              measure_one_at_time = FALSE,
                              minimize_only = FALSE,
                              data = NULL,
                              max_steps = 1000,
                              return_step_results = FALSE,
                              do_not_remove = c("(Intercept)")) {
  log_ <- list()
  cur_step <- 0
  seen_states <- list(names(stats::coef(model)))
  do_backward <- direction %in% c("backward", "both")
  do_forward <- direction %in% c("forward", "both")

  threshold_fn <- if (is.function(threshold)) {
    threshold
  } else {
    function(model_) threshold
  }

  if (is.null(data)) data <- stats::model.frame(model)
  while (cur_step < max_steps) {
    cur_step <- cur_step + 1

    cur_threshold <- threshold_fn(model)

    # 1. Backward removal
    if (do_backward && length(stats::coef(model)) > 1) {
      values <- backward_values(
        model = model,
        measure_fn = measure_fn,
        measure_one_at_time = measure_one_at_time,
        data = data,
        do_not_remove = do_not_remove,
        seen_states = seen_states
      )

      updated_model <- update_model_remove(
        model = model,
        values = values,
        threshold = cur_threshold,
        data = data,
        minimize_only = minimize_only
      )

      if (!is.null(updated_model)) {
        model <- updated_model$fit
        seen_states <- append(seen_states, list(names(stats::coef(model))))
        log_[[cur_step]] <- list(
          step = cur_step, operation = "removal",
          value = updated_model$removed_var
        )
        next
      }
    }

    # 2. Forward selection
    if (do_forward) {
      values <- forward_values(
        model = model,
        measure_fn = measure_fn,
        measure_one_at_time = measure_one_at_time,
        data = data,
        addable_coefs = addable_coefs,
        seen_states = seen_states
      )

      updated_model <- update_model_add(model, values, cur_threshold, data)
      if (!is.null(updated_model)) {
        model <- updated_model$fit
        seen_states <- append(seen_states, list(names(stats::coef(model))))
        log_[[cur_step]] <- list(
          step = cur_step, operation = "addition",
          value = updated_model$added_var
        )
        next
      }
    }

    # If cant remove or add a variable stop the selection
    break
  }

  if (return_step_results) {
    return(list(fit = model, log = log_))
  }

  return(model)
}
