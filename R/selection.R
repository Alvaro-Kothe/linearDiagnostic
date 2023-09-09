remove_variable <- function(x) {

  INTERCEPT = "(Intercept)"
  remove_intercept <- INTERCEPT %in% x

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
#' argument the selection result, and second argument a named vector with
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
  removed_values <- double()
  removed_names <- character()
  cur_step <- 0

  while (cur_step < max_steps) {
    cur_step <- cur_step + 1
    values <- measure_fn(model)
    explore_remove <- setdiff(names(values), do_not_remove)
    values <- values[explore_remove]
    to_remove <- which.max(values)

    if (values[to_remove] < threshold) {
      break
    }

    removed_values[cur_step] <- values[to_remove]
    removed_names[cur_step] <- removed_var <- names(to_remove)
    next_formula <- remove_variable(removed_var)

    model <- stats::update(model, formula. = next_formula, data = data)
  }
  removed <- stats::setNames(removed_values, removed_names)

  if (return_step_results) return(list(fit = model, removed = removed))

  return(model)
}
