collapse_vec <- function(x) {
  base::sprintf("c(%s)", paste0(x, collapse = ", "))
}

change_reponse_formula <- function(x) {
  formula_string <- if (is.matrix(x)) {
    columns_collapsed <- apply(x, 2, collapse_vec, simplify = FALSE)
    base::sprintf("cbind(%s) ~ .", paste(columns_collapsed, collapse = ", "))
  } else {
    base::sprintf("%s ~ .", collapse_vec(x))
  }
  stats::as.formula(formula_string)
}

#' Refit a model
#'
#' Refit a model with a new response.
#'
#' This function uses `new_response` to refit `object` replacing its old response variable.
#' If the class is `merMod` it uses [lme4::refit()], otherwise uses [stats::update()].
#'
#' The default method tries to update the model response using it's [stats::model.frame()],
#' if it errors it tries to update the model by inserting the `new_response`
#' directly into the object formula.
#'
#' @param object A model.
#' @param new_response the new response, may be a vector or a matrix.
#' @param ... other arguments passed to `refit` or `update`.
#' @return A model with same class as `object`.
#' @seealso [stats::update()]
#' @export
get_refit <- function(object, new_response, ...) {
  UseMethod("get_refit", object)
}

#' @rdname get_refit
#' @export
get_refit.default <- function(object, new_response, ...) {
  if (!is.vector(new_response) && !is.matrix(new_response)) {
    stop("`new_response` should be either a vector or matrix")
  }
  if (as.character(stats::formula(object)[[2]])[[1]] == "c") {
    # The model was fitted using c(vec) ~ X
    return(update_using_formula(object, new_response, ...))
  }
  tryCatch(
    update_using_model_frame(object, new_response, ...),
    error = function(e) update_using_formula(object, new_response, ...)
  )
}

#' Update the object with new response using only model frame.
#'
#' @inheritParams get_refit
update_using_model_frame <- function(object, new_response, ...) {
  model_frame <- stats::model.frame(object)
  if (is.vector(new_response)) {
    model_frame[, 1] <- new_response
  } else if (is.matrix(new_response)) {
    model_frame[[1]] <- NULL
    model_frame <- cbind(new_response, model_frame)
  }
  stats::update(object, data = model_frame, ...)
}

#' Update the object using only the formula
#'
#' Create a new formula with `new_response` values directly on the formula left hand side.
#'
#' This function is designed to be used as a fallback for [update_using_model_frame()]
#' as it's prone to run with memory issues.
#'
#' @inheritParams get_refit
update_using_formula <- function(object, new_response, ...) {
  new_formula <- change_reponse_formula(new_response)
  stats::update(object, formula. = new_formula, ...)
}

#' @rdname get_refit
#' @export
get_refit.merMod <- function(object, new_response, ...) {
  lme4::refit(object, new_response, ...)
}
