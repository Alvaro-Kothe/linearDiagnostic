collapse_vec <- function(x) {
  base::sprintf("c(%s)", paste0(x, collapse = ", "))
}

change_reponse_formula <- function(x) {
  formula_string <- if ("matrix" %in% class(x)) {
    columns_collapsed <- apply(x, 2, collapse_vec, simplify = FALSE)
    base::sprintf("cbind(%s) ~ .", paste(columns_collapsed, collapse = ", "))
  } else {
    base::sprintf("%s ~ .", collapse_vec(x))
  }
  stats::as.formula(formula_string)
}
