test_that("envelope() observed and expected concides with qqnorm", {
  fit <- lm(mpg ~ cyl, data = mtcars)
  residual_fn <- stats::rstudent
  env_meas <- envelope(fit, residual_fn = residual_fn, n_sim = 5, plot.it = FALSE)
  residual <- residual_fn(fit)
  expect_mapequal(env_meas$observed, residual)
  expect_mapequal(env_meas$observed, stats::qqnorm(residual, plot.it = FALSE)$y)
  expect_setequal(env_meas$expected, stats::qqnorm(residual, plot.it = FALSE)$x)
})

test_that("envelope() acceptable coverage for correct model", {
  set.seed(1)
  x <- matrix(rnorm(1000), ncol = 2)
  beta <- c(1, 1)
  y <- rnorm(nrow(x), mean = x %*% beta, sd = 0.01)
  fit <- lm(y ~ x + 0)
  env_meas <- envelope(fit, n_sim = 10, plot.it = FALSE)
  inside_band <- mean(env_meas$lower < env_meas$observed & env_meas$observed < env_meas$upper)
  expect_gte(inside_band, 0.87)
  expect_equal(1 - inside_band, mean(env_meas$outside))
})

test_that("envelope() detects incorrect fit", {
  set.seed(1)
  x <- matrix(rnorm(1000), ncol = 2)
  beta <- c(1, 1)
  y <- rpois(nrow(x), lambda = exp(x %*% beta))
  fit <- lm(y ~ x + 0)
  env_meas <- envelope(fit, n_sim = 10, plot.it = FALSE)
  outside_band <- mean(env_meas$outside)
  expect_gt(outside_band, 0.5)
})

test_that("plot envelope runs without errors", {
  expect_no_error({
    env_meas <- envelope(lm(c(1, 5) ~ 1), n_sim = 2, plot.it = FALSE)
    plot(env_meas)
  })
})

test_that("envelope() is compatible with models using cbind", {
  m <- c(1, 4, 10, 30)
  y <- c(0, 2, 5, 15)
  fit <- glm(cbind(y, m - y) ~ 1, family = binomial())
  expect_no_error(envelope(fit, n_sim = 2, plot.it = FALSE))
})

test_that("envelope works with lme4::lmer", {
  data("sleepstudy", package = "lme4")
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  expect_no_error(envelope(fit, n_sim = 2, residual_fn = residuals, plot.it = FALSE))
})

test_that("envelope works with lme4::glmer", {
  data("cbpp", package = "lme4")
  fit <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp, family = binomial
  )
  expect_no_error(envelope(fit, n_sim = 2, residual_fn = residuals, plot.it = FALSE))
})

test_that("envelope works with glmmTMB", {
  data("Salamanders", package = "glmmTMB")
  m1 <- glmmTMB::glmmTMB(count ~ mined + (1 | site),
    zi = ~mined,
    family = poisson, data = Salamanders
  )
  expect_no_error(envelope(m1, n_sim = 2, residual_fn = residuals))
  m2 <- glmmTMB::glmmTMB(count ~ spp + mined + (1 | site),
    zi = ~ spp + mined,
    family = glmmTMB::nbinom2, data = Salamanders
  )
  expect_no_error(envelope(m2, n_sim = 2, residual_fn = residuals))
})
