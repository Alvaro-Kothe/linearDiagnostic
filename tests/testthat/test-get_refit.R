test_that("response is replaced when using c() as response", {
  fit <- lm(c(2, 4, 5) ~ c(3, 5, 7))
  model_refit <- get_refit(fit, c(1, 1, 1))
  expect_equal(stats::model.frame(model_refit)[[1]], c(1, 1, 1))
})

test_that("response is replaced when using data", {
  df <- data.frame(y = c(2, 4, 5), x = c(1, 2, 3))
  fit <- lm(y ~ x, data = df)
  model_refit <- get_refit(fit, c(1, 1, 1))
  expect_equal(stats::model.frame(model_refit)[[1]], c(1, 1, 1))
})
