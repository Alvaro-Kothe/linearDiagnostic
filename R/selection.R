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

get_p_value <- function(model) {
  return(summary(model)[["coefficients"]][, 4])
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
                               measure_fn = get_p_value,
                               data = stats::model.frame(model),
                               max_steps = 1000,
                               return_step_results = FALSE,
                               do_not_remove = c("(Intercept)")) {
  log_ <- list()
  cur_step <- 0

  while (cur_step < max_steps) {
    cur_step <- cur_step + 1
    values <- measure_fn(model)
    explore_remove <- setdiff(names(values), do_not_remove)
    values <- values[explore_remove]
    to_remove <- which.max(values)

    if (length(values) == 0L || values[to_remove] < threshold) {
      break
    }

    removed_var <- names(to_remove)
    log_[[cur_step]] <- list(
      step = cur_step, operation = "removal",
      variable = removed_var, value = values[to_remove]
    )
    next_formula <- remove_variable(removed_var)

    model <- stats::update(model, formula. = next_formula, data = data)
  }

  if (return_step_results) {
    return(list(fit = model, log = log_))
  }

  return(model)
}


#' Select covariates bidirectional
#'
#' @param model A model with [stats::update()] method.
#' @param threshold Value threshold to remove variable. where the variable is
#' removed if `measure_fn(model) > threshold`, added if `measure_fn(model) <= threshold`
#' @param addable_coefs Coefficients to try to add during forward step.
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
#' bidirectional_selection(model)
bidirectional_selection <- function(model, threshold = .15,
                                    addable_coefs = NULL,
                                    measure_fn = get_p_value,
                                    data = stats::model.frame(model),
                                    max_steps = 1000,
                                    return_step_results = FALSE,
                                    do_not_remove = c("(Intercept)")) {
  log_ <- list()
  cur_step <- 0
  seen_states <- list(names(coef(model)))
  if (is.null(addable_coefs)) addable_coefs <- names(coef(model))

  if (is.null(data)) data <- model.frame(model)
  while (cur_step < max_steps) {
    cur_step <- cur_step + 1

    # 1. Backward removal
    values <- measure_fn(model)

    possible_next_states <- lapply(names(values), function(x) {
      names(coef(model))[names(coef(model)) != x]
    })
    next_state_seen <- sapply(
      possible_next_states,
      function(x) any(sapply(seen_states, setequal, x))
    )
    values <- values[!next_state_seen]

    explore_remove <- setdiff(names(values), do_not_remove)
    values <- values[explore_remove]

    to_remove <- which.max(values)

    if (length(values) > 0 && values[to_remove] > threshold) {
      removed_var <- names(to_remove)
      next_formula <- remove_variable(removed_var)

      model <- update(model, formula. = next_formula, data = data)

      seen_states <- append(seen_states, list(names(coef(model))))
      log_[[cur_step]] <- list(
        step = cur_step, operation = "removal",
        variable = removed_var, value = values[to_remove]
      )
      next
    }

    # 2. Forward selection
    add_candidates <- setdiff(addable_coefs, names(coef(model)))
    possible_next_states <- lapply(add_candidates, c, names(coef(model)))
    next_state_seen <- sapply(
      possible_next_states,
      function(x) any(sapply(seen_states, setequal, x))
    )
    add_candidates <- add_candidates[!next_state_seen]

    values <- sapply(add_candidates, function(add_cand) {
      next_formula <- add_variable(add_cand)

      next_fit <- update(model, formula. = next_formula, data = data)

      measure_fn(next_fit)[[add_cand]]
    })

    to_add <- which.min(values)

    if (length(values) > 0 && values[to_add] <= threshold) {
      added_var <- names(to_add)
      next_formula <- add_variable(added_var)

      model <- update(model, formula. = next_formula, data = data)

      seen_states <- append(seen_states, list(names(coef(model))))
      log_[[cur_step]] <- list(
        step = cur_step, operation = "addition",
        variable = added_var, value = values[to_add]
      )
      next
    }

    # If cant remove or add a variable stop the selection
    break
  }

  if (return_step_results) {
    return(list(fit = model, log = log_))
  }

  return(model)
}
