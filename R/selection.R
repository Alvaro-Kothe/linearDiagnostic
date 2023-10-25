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
#' @param model A model with [stats::update()] method.
#' @param threshold Value threshold to remove variable. where the variable is
#' removed if `measure_fn(model) >= threshold`.
#' @param measure_fn Function with model as argument and returns values to be used by
#' `threshold`. Default is p-value reported by [summary()].
#' @param data Data to be used for model refit.
#' @param max_steps Maximum number of steps for selection.
#' @param return_step_results Boolean that if `TRUE` returns a list with the first
#' argument the selection result, and second argument a list with step, operation,
#' variable and its value.
#' `measure_fn` measures and coefficient removed.
#' @param do_not_remove String vector of variables to not be removed in the selection.
#'
#' @return the selected model if `return_step_results` is `False`.
#' @export
#'
#' @examples
#' model <- lm(mpg ~ ., data = mtcars)
#' backward_selection(model)
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
#' @param model A model with [stats::update()] method.
#' @param threshold Value threshold to remove variable. where the variable is
#' removed if `measure_fn(model) > threshold`, added if `measure_fn(model) <= threshold`
#' @param addable_coefs Coefficients to try to add during forward step.
#' @param measure_fn Function with model as argument and returns values to be used by
#' `threshold`. It can also compare two models, where during forward step
#' it calls `measure_fn(candidate_model, current_selected_model)` and
#' during backward step it calls `measure_fn(current_selected_model, candidate_model)`.
#' Default is p-value reported by [summary()].
#' @param data Data to be used for model refit.
#' @param max_steps Maximum number of steps for selection.
#' @param return_step_results Boolean that if `TRUE` returns a list with the first
#' argument the selection result, and second argument a list with step, operation,
#' variable and its value.
#' `measure_fn` measures and coefficient removed.
#' @param do_not_remove String vector of variables to not be removed in the selection.
#'
#' @return the selected model if `return_step_results` is `False`.
#' @export
#'
#' @examples
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
                            data,
                            do_not_remove, seen_states) {
  nargs_measure_fn <- length(formals(measure_fn))
  cur_coefs <- names(stats::coef(model))
  candidates_remove <- stats::setNames(nm = setdiff(cur_coefs, do_not_remove))


  possible_next_states <- lapply(candidates_remove, function(x) {
    names(stats::coef(model))[names(stats::coef(model)) != x]
  })

  next_state_seen <- sapply(
    possible_next_states,
    function(x) any(sapply(seen_states, setequal, x))
  )

  next_states_lookup <- names(next_state_seen)[!next_state_seen]

  if (nargs_measure_fn == 1L) {
    values <- measure_fn(model)[next_states_lookup]
  } else if (nargs_measure_fn == 2L) {
    values <- sapply(next_states_lookup, function(remove_candidate) {
      next_formula <- remove_variable(remove_candidate)

      next_fit <- stats::update(model, formula. = next_formula, data = data)

      return(measure_fn(model, next_fit))
    })
  }

  return(values)
}

forward_values <- function(model, measure_fn, data, addable_coefs, seen_states) {
  nargs_measure_fn <- length(formals(measure_fn))
  add_candidates <- setdiff(addable_coefs, names(stats::coef(model)))
  possible_next_states <- lapply(add_candidates, c, names(stats::coef(model)))
  next_state_seen <- sapply(
    possible_next_states,
    function(x) any(sapply(seen_states, setequal, x))
  )
  add_candidates <- add_candidates[!next_state_seen]

  values <- sapply(add_candidates, function(add_cand) {
    next_formula <- add_variable(add_cand)

    next_fit <- stats::update(model, formula. = next_formula, data = data)

    if (nargs_measure_fn == 2L) {
      return(measure_fn(next_fit, model))
    }

    return(measure_fn(next_fit)[[add_cand]])
  })

  return(values)
}

update_model_remove <- function(model, values, threshold, data) {
  to_remove <- which.max(values)
  if (length(to_remove) == 0L) {
    return(NULL)
  }

  value <- values[to_remove]

  if (value <= threshold) {
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
#' @param model A model with [stats::update()] method.
#' @param threshold Value threshold to remove variable. where the variable is
#' removed if `measure_fn(model) > threshold`, added if `measure_fn(model) <= threshold`
#' @param direction Selection method, being both(forward and backward), forward and backward.
#' default is "both"
#' @param addable_coefs Coefficients to try to add during forward step.
#' @param measure_fn Function with model as argument and returns values to be used by
#' `threshold`. It can also compare two models, where during forward step
#' it calls `measure_fn(candidate_model, current_selected_model)` and
#' during backward step it calls `measure_fn(current_selected_model, candidate_model)`.
#' @param data Data to be used for model refit.
#' @param max_steps Maximum number of steps for selection.
#' @param return_step_results Boolean that if `TRUE` returns a list with the first
#' argument the selection result, and second argument a list with step, operation,
#' variable and its value.
#' `measure_fn` measures and coefficient removed.
#' @param do_not_remove String vector of variables to not be removed in the selection.
#'
#' @return the selected model if `return_step_results` is `False`.
#' @export
#'
#' @examples
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
select_covariates <- function(model,
                              threshold = .15,
                              direction = "both",
                              addable_coefs = names(stats::coef(model)),
                              measure_fn = function(x) summary(x)[["coefficients"]][, 4],
                              data = NULL,
                              max_steps = 1000,
                              return_step_results = FALSE,
                              do_not_remove = c("(Intercept)")) {
  log_ <- list()
  cur_step <- 0
  seen_states <- list(names(stats::coef(model)))
  do_backward <- direction %in% c("backward", "both")
  do_forward <- direction %in% c("forward", "both")

  if (is.null(data)) data <- stats::model.frame(model)
  while (cur_step < max_steps) {
    cur_step <- cur_step + 1

    # 1. Backward removal
    if (do_backward && length(stats::coef(model)) > 1) {
      values <- backward_values(model, measure_fn, data, do_not_remove, seen_states)
      updated_model <- update_model_remove(model, values, threshold, data)
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
      values <- forward_values(model, measure_fn, data, addable_coefs, seen_states)

      updated_model <- update_model_add(model, values, threshold, data)
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
