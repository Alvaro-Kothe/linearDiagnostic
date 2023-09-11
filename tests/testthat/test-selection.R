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
  selection <- bidirectional_selection(fit)
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
