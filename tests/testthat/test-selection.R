test_that("backward_selection() without covariates works", {
  fit <- lm(mpg ~ 1, data = mtcars)
  selection <- backward_selection(fit)
  expect_equal(names(coef(selection)), "(Intercept)")
})

test_that("backward_selection() remove covariates", {
  fit <- lm(mpg ~ ., data = mtcars)
  selection <- backward_selection(fit)
  expect_lt(length(coef(selection)), length(coef(fit)))
})

test_that("backward_selection() removes Intercept", {
  set.seed(1)
  x <- c(1, 2, 5, 2, 3, 4)
  y <- 2 * x + runif(length(x), .01)
  fit <- lm(y ~ x)
  selection <- backward_selection(fit, do_not_remove = NULL)
  expect_false("(Intercept)" %in% names(coef(selection)))
})

test_that("backward_selection() keeps covariates from `do_not_remove`", {
  fit <- lm(mpg ~ ., data = mtcars)
  keep <- c("(Intercept)", "cyl", "am")
  selection <- backward_selection(fit,
    data = mtcars,
    do_not_remove = keep
  )
  expect_in(keep, names(coef(selection)))
})

test_that("bidirectional_selection() keeps covariates from `do_not_remove`", {
  fit <- lm(mpg ~ ., data = mtcars)
  keep <- c("(Intercept)", "cyl", "am")
  selection <- bidirectional_selection(fit,
    data = mtcars,
    do_not_remove = keep
  )
  expect_in(keep, names(coef(selection)))
})

test_that("bidirectional_selection() add covariates", {
  fit <- lm(mpg ~ 1, data = mtcars)
  selection <- bidirectional_selection(fit,
    addable_coefs = setdiff(colnames(mtcars), "mpg"),
    data = mtcars
  )
  expect_gt(length(coef(selection)), 1)
})

test_that("bidirectional_selection() remove covariates", {
  fit <- lm(mpg ~ ., data = mtcars)
  selection <- bidirectional_selection(fit, data = mtcars)
  expect_lt(length(coef(selection)), length(coef(fit)))
})

test_that("bidirectional_selection() removes Intercept", {
  set.seed(1)
  x <- c(1, 2, 5, 2, 3, 4)
  y <- 2 * x + runif(length(x), .01)
  fit <- lm(y ~ x)
  selection <- bidirectional_selection(fit, do_not_remove = NULL)
  expect_false("(Intercept)" %in% names(coef(selection)))
})

test_that("Likelihood ratio test works for selection", {
  fit <- lm(mpg ~ ., data = mtcars)

  lrt <- function(model1, model2) {
    lrt_stat <- 2 * (logLik(model1)[1L] - logLik(model2)[1L])
    return(1 - pchisq(lrt_stat, 1))
  }

  selection <- select_covariates(fit, measure_fn = lrt)

  expect_false(setequal(
    names(coef(fit)),
    names(coef(selection))
  ))
})

test_that("Likelihood ratio test works for forward selection", {
  fit <- lm(mpg ~ 1, data = mtcars)

  lrt <- function(model1, model2) {
    lrt_stat <- 2 * (logLik(model1)[1L] - logLik(model2)[1L])
    return(1 - pchisq(lrt_stat, 1))
  }

  selection <- select_covariates(fit,
    measure_fn = lrt,
    direction = "forward",
    addable_coefs = setdiff(colnames(mtcars), "mpg"),
    data = mtcars
  )

  expect_gt(length(coef(selection)), 1)
})

test_that("Selection works with AIC", {
  fit <- lm(mpg ~ ., data = mtcars)

  selection <- select_covariates(fit,
    measure_fn = function(model) stats::extractAIC(model)[2L],
    threshold = function(model) stats::extractAIC(model)[2L],
    measure_one_at_time = TRUE,
    direction = "both",
    minimize_only = TRUE,
    data = mtcars
  )

  selection_step_aic <- step(fit, trace = 0)

  expect_equal(coef(selection), coef(selection_step_aic))
})

test_that("Selection works with AIC start from intercept", {
  fit <- lm(mpg ~ 1, data = mtcars)
  upper_scope <- ~ cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb

  selection <- select_covariates(fit,
    measure_fn = function(model) stats::extractAIC(model)[2L],
    threshold = function(model) stats::extractAIC(model)[2L],
    addable_coefs = setdiff(colnames(mtcars), "mpg"),
    measure_one_at_time = TRUE,
    direction = "both",
    minimize_only = TRUE,
    data = mtcars
  )

  selection_step_aic <- step(
    fit,
    scope = list(upper = upper_scope, lower = ~1),
    trace = 0
  )

  expect_equal(coef(selection), coef(selection_step_aic))
})

test_that("return_step_results works", {
  set.seed(1)
  x <- c(1, 2, 5, 2, 3, 4)
  y <- 2 * x + runif(length(x), .01)
  fit <- lm(y ~ x)
  expect_no_error(bidirectional_selection(fit, return_step_results = TRUE)$log)
  expect_no_error(backward_selection(fit, return_step_results = TRUE)$log)
})

test_that("update_model_*() returns expected output", {
  fit <- lm(mpg ~ cyl + disp, data = mtcars)
  expect_null(update_model_remove(fit, double(), .1, data = mtcars))
  expect_null(update_model_remove(fit, c(cyl = .05), .1, data = mtcars))
  expect_null(update_model_add(fit, double(), .1, data = mtcars))
  expect_null(update_model_add(fit, c(hp = .5), .1, data = mtcars))

  updated_remove <- update_model_remove(fit, c(cyl = .2), .1, data = mtcars)
  expect_true(setdiff(
    names(coef(fit)),
    names(coef(updated_remove$fit))
  ) == c("cyl"))

  expect_identical(updated_remove$removed_var, c(cyl = .2))

  updated_add <- update_model_add(fit, c(hp = .1), .2, data = mtcars)
  expect_true(setdiff(
    names(coef(updated_add$fit)),
    names(coef(fit))
  ) == c("hp"))

  expect_identical(updated_add$added_var, c(hp = .1))
})

test_that("selection errors with measure_fn with 0 or 3+ arguments", {
  fit <- lm(mpg ~ ., data = mtcars)
  expect_error(select_covariates(fit, measure_fn = function() {
    42
  }))
  expect_error(select_covariates(fit, measure_one_at_time = TRUE, measure_fn = function(a, b, c) {
    42
  }))
})
