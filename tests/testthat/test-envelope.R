test_that("envelope_measures() observed and expected concides with qqnorm", {
  fit <- lm(mpg ~ cyl, data = mtcars)
  residual_fn <- stats::rstudent
  env_meas <- envelope_measures(fit, residual_fn = residual_fn)
  residual <- residual_fn(fit)
  expect_mapequal(env_meas$observed, residual)
  expect_mapequal(env_meas$observed, stats::qqnorm(residual, plot.it = FALSE)$y)
  expect_setequal(env_meas$expected, stats::qqnorm(residual, plot.it = FALSE)$x)
})

test_that("envelope_measures() acceptable coverage for correct model", {
  x <- matrix(rnorm(1000), ncol = 2)
  beta <- c(1, 1)
  y <- rnorm(nrow(x), mean = x %*% beta, sd = .01)
  fit <- lm(y ~ x + 0)
  env_meas <- envelope_measures(fit)
  inside_band <- mean(env_meas$lower < env_meas$observed &
    env_meas$observed < env_meas$upper)
  expect_gte(inside_band, .90)
  expect_equal(1 - inside_band, mean(env_meas$outside))
})

test_that("envelope_measures() detects incorrect fit", {
  x <- matrix(rnorm(1000), ncol = 2)
  beta <- c(1, 1)
  y <- rpois(nrow(x), lambda = exp(x %*% beta))
  fit <- lm(y ~ x + 0)
  env_meas <- envelope_measures(fit)
  outside_band <- mean(env_meas$outside)
  expect_gt(outside_band, .5)
})
