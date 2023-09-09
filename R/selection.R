remove_variable <- function(x) {
  rhs <- paste0(x, collapse = "-")
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
