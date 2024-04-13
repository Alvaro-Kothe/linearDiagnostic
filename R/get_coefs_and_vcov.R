#' Get fixed effects
#'
#' Extracts the fixed effects coefficients from a model.
#'
#' By default it calls [stats::coef()]. If the model class is `merMod` calls
#' `fixef`.
#'
#' @param object A model object for which fixed effects are to be retrieved.
#' @return A numeric vector of fixed effects coefficients.
#'
#' @export
get_fixef <- function(object) {
  UseMethod("get_fixef")
}

#' @rdname get_fixef
#' @export
get_fixef.default <- function(object) {
  stats::coef(object)
}

#' @rdname get_fixef
#' @export
get_fixef.merMod <- function(object) {
  lme4::fixef(object)
}

#' Get covariance matrix
#'
#' Retrieves the covariance matrix for the fixed effects.
#'
#' By default it calls [stats::vcov()] to retrieve the covariances.
#'
#' @param object A model object for which the covariance matrix is to be retrieved.
#' @return A covariance matrix of the model object.
#' @export
get_vcov <- function(object) {
  UseMethod("get_vcov")
}

#' @rdname get_vcov
#' @export
get_vcov.default <- function(object) {
  as.matrix(stats::vcov(object))
}
