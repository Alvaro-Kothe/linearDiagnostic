% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection.R
\name{select_covariates}
\alias{select_covariates}
\title{Select covariates}
\usage{
select_covariates(
  model,
  threshold = 0.15,
  direction = "both",
  addable_coefs = names(stats::coef(model)),
  measure_fn = function(x) summary(x)[["coefficients"]][, 4],
  measure_one_at_time = FALSE,
  minimize_only = FALSE,
  data = NULL,
  max_steps = 1000,
  return_step_results = FALSE,
  do_not_remove = c("(Intercept)")
)
}
\arguments{
\item{model}{A model with [stats::update()] method.}

\item{threshold}{Value threshold to remove variable. where the variable is
removed if `measure_fn(model) > threshold`, added if `measure_fn(model) <= threshold`.
It can also be a function that evaluates on the current fitted model.}

\item{direction}{Selection method, being both(forward and backward), forward and backward.
default is "both"}

\item{addable_coefs}{Coefficients to try to add during forward step.}

\item{measure_fn}{Function with model as argument and returns values to be used by
`threshold`. It can also compare two models, where during forward step
it calls `measure_fn(candidate_model, current_selected_model)` and
during backward step it calls `measure_fn(current_selected_model, candidate_model)`.}

\item{measure_one_at_time}{Boolean indicating to apply `measure_fn` in one model at a time.
Set this option to `TRUE` if `measure_fn` returns an atomic value, for example if
`measure_fn` is `AIC`.}

\item{minimize_only}{Boolean indicating that during backward model update
it should minimize the `measure_fn` instead of maximize it.}

\item{data}{Data to be used for model refit.}

\item{max_steps}{Maximum number of steps for selection.}

\item{return_step_results}{Boolean that if `TRUE` returns a list with the first
argument the selection result, and second argument a list with step, operation,
variable and its value.
`measure_fn` measures and coefficient removed.}

\item{do_not_remove}{String vector of variables to not be removed in the selection.}
}
\value{
the selected model if `return_step_results` is `False`.
}
\description{
Select covariates
}
\examples{
model <- lm(mpg ~ ., data = mtcars)
select_covariates(model)

## measure_fn with two parameters

lrt <- function(model1, model2) {
  lrt_stat <- 2 * (logLik(model1)[1L] - logLik(model2)[1L])
  return(1 - pchisq(lrt_stat, 1))
}

select_covariates(model, measure_fn = lrt)

## AICc selection

AICc <- function(model) {
  loglike <- logLik(model)
  df <- attr(loglike, "df")
  nobs <- attr(loglike, "nobs")
  aic <- -2 * as.numeric(loglike) + 2 * df

  aicc <- aic + (2 * (df^2) + 2 * df) / (nobs - df - 1)

  return(aicc)
}

selection <- select_covariates(model,
  measure_fn = AICc,
  threshold = AICc,
  measure_one_at_time = TRUE,
  minimize_only = TRUE,
  direction = "both",
  data = mtcars
)
}
