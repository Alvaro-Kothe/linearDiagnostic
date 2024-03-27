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
  set.seed(1)
  x <- matrix(rnorm(1000), ncol = 2)
  beta <- c(1, 1)
  y <- rnorm(nrow(x), mean = x %*% beta, sd = .01)
  fit <- lm(y ~ x + 0)
  env_meas <- envelope_measures(fit)
  inside_band <- mean(env_meas$lower < env_meas$observed & env_meas$observed < env_meas$upper)
  expect_gte(inside_band, .87)
  expect_equal(1 - inside_band, mean(env_meas$outside))
})

test_that("envelope_measures() detects incorrect fit", {
  set.seed(1)
  x <- matrix(rnorm(1000), ncol = 2)
  beta <- c(1, 1)
  y <- rpois(nrow(x), lambda = exp(x %*% beta))
  fit <- lm(y ~ x + 0)
  env_meas <- envelope_measures(fit)
  outside_band <- mean(env_meas$outside)
  expect_gt(outside_band, .5)
})

test_that("plot_envelope() runs without errors", {
  expect_no_error(plot_envelope(lm(c(1, 5) ~ 1)))
})

test_that("plot_envelope() is compatible with models using cbind", {
  m <- c(1, 4, 10, 30)
  y <- c(0, 2, 5, 15)
  fit <- glm(cbind(y, m - y) ~ 1, family = binomial())
  expect_no_error(plot_envelope(fit))
})
