simple_lm_fit <- function() {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  lm(y ~ x)
}

test_that("set.seed generates same stuff", {
  fit <- simple_lm_fit()

  set.seed(1)
  simulation_result1 <- get_p_values(fit, n_sim = 5)

  set.seed(1)
  simulation_result2 <- get_p_values(fit, n_sim = 5)

  expect_identical(simulation_result1, simulation_result2)
})

test_that("test_coefficients generate different results", {
  fit <- simple_lm_fit()
  set.seed(1)
  p_values_regular <- get_p_values(fit, n_sim = 10)

  test_coefs <- c(1000000, 1000000)
  set.seed(1)
  p_values_test_coefs <- get_p_values(fit, n_sim = 10, test_coefficients = test_coefs)

  # Test if all elements are close to 0, as test_coefs were very high
  expect_true(all(p_values_test_coefs$pvalues_joint < 1e-8))
  expect_true(all(p_values_test_coefs$pvalues_matrix < 1e-8))
  expect_false(identical(p_values_regular, p_values_test_coefs))
})

test_that("rownames matrix is the same as model's coefficients names", {
  df <- data.frame("foo" = c(1, 2, 3, 7), "bar" = c(2, 3, 1, 9), "y" = c(3, 5, 1, 2))
  fit <- lm(y ~ foo + bar, data = df)
  p_values <- get_p_values(fit, n_sim = 2)
  expect_named(p_values$pvalues_matrix[, 1], c("(Intercept)", "foo", "bar"))
})

test_that("generator yield different results", {
  fit <- simple_lm_fit()

  suppressWarnings({
    set.seed(1)
    simulation_result1 <- get_p_values(fit, n_sim = 5)
    set.seed(1)
    simulation_result2 <- get_p_values(fit,
      n_sim = 5,
      generator = function(model) rgamma(stats::nobs(model), abs(stats::fitted.values(model)))
    )
  })
  expect_false(identical(simulation_result1, simulation_result2))
})

test_that("Same response yield same result", {
  fit <- simple_lm_fit()
  set.seed(1)
  responses <- simulate(fit, 5)
  sim1 <- get_p_values(fit, n_sim = 5, responses = responses)
  sim2 <- get_p_values(fit, n_sim = 5, responses = responses)
  expect_identical(sim1, sim2)
})

test_that("p_values are equal with deterministic generator", {
  fit <- simple_lm_fit()
  generator <- function(...) 1:4
  suppressWarnings(
    p_values <- get_p_values(fit, n_sim = 5, generator = generator)
  )
  expect_setequal(p_values$pvalues_joint, p_values$pvalues_joint[1])
  expect_setequal(p_values$pvalues_matrix[1, ], p_values$pvalues_matrix[1, 1])
  expect_setequal(p_values$pvalues_matrix[2, ], p_values$pvalues_matrix[2, 1])
})

test_that("get_p_values throw warning with singular matrix", {
  params <- c(1, .2, .5, -.2, -.5, 0, 0.1, 0.01, 3, 5)
  set.seed(1)
  model_matrix <- matrix(rnorm(13 * 10),
    nrow = 13,
    ncol = 10
  )
  eta <- model_matrix %*% params
  mu <- exp(eta)
  y <- rpois(13, mu)

  suppressWarnings({
    fit <- glm(y ~ model_matrix + 1, family = poisson())

    expect_warning(
      get_p_values(fit, n_sim = 10),
      "Couldn't inverse vcov from \\d+ simulations and used `MASS::ginv` instead"
    )
  })
})

test_that("get_p_values works with poisson with offset", {
  n <- 5
  offset_ <- rpois(n, 10000)
  x <- runif(n)
  y <- rpois(n, exp(1 + 2 * x))

  suppressWarnings({
    fit <- glm(y ~ x + offset(offset_), family = poisson())
    p_values <- get_p_values(fit, n_sim = 5)
  })

  expect_length(p_values$simulation_fixef[[1]], 2)
})

test_that("get_p_values convergence for lm is NA", {
  fit <- simple_lm_fit()
  p_values <- get_p_values(fit, n_sim = 2)
  expect_equal(p_values$converged, c(NA, NA))
})

test_that("get_p_values convergence for glm is logical", {
  fit <- glm(c(1, 3, 5) ~ c(1, 2, 3), family = poisson())
  set.seed(1)
  p_values <- get_p_values(fit, n_sim = 2)
  expect_equal(p_values$converged, c(TRUE, TRUE))
})

test_that("Can compute p_values from merMod class", {
  data("sleepstudy", package = "lme4")
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  set.seed(1)
  expect_no_error(p_values <- get_p_values(fit, n_sim = 2))
  expect_equal(p_values$converged, c(TRUE, TRUE))
})

test_that("get_p_values() with custom method works", {
  foo <- function(formula, data) {
    model_frame <- model.frame(formula, data = data)
    y <- model.response(model_frame, type = "numeric")
    x <- model.matrix(formula, data)
    xtx <- crossprod(x)
    xtxinv <- solve(xtx)

    out <- list()
    out$coef <- c(xtxinv %*% t(x) %*% y)
    out$fitted.values <- c(x %*% out$coef)
    out$sigma <- sum((y - out$fitted.values)^2) / (length(y) - length(out$coef))
    out$x <- x
    out$y <- y
    out$data <- data
    out$model_frame <- model_frame
    out$vcov <- out$sigma * xtxinv
    out$formula <- formula

    class(out) <- "foo"

    return(out)
  }

  # nolint start: object_name_linter
  nobs.foo <- function(object, ...) {
    length(object$y)
  }
  assign("nobs.foo", nobs.foo, envir = .GlobalEnv)

  simulate.foo <- function(object, nsim = 1, seed = NULL, ...) {
    mu <- object$fitted.values

    replicate(nsim, rep(NA, length(mu)), simplify = FALSE)
  }
  assign("simulate.foo", simulate.foo, envir = .GlobalEnv)

  update.foo <- function(object, ...) {
    object
  }
  assign("update.foo", update.foo, envir = .GlobalEnv)

  vcov.foo <- function(object) {
    matrix(NA, nrow = 2, ncol = 2)
  }
  assign("vcov.foo", vcov.foo, envir = .GlobalEnv)

  coef.foo <- function(object) {
    rep_len(NA, 2)
  }
  assign("coef.foo", coef.foo, envir = .GlobalEnv)

  formula.foo <- function(object) {
    object$formula
  }
  assign("formula.foo", formula.foo, envir = .GlobalEnv)
  # nolint end


  fit <- foo(mpg ~ cyl, mtcars)

  suppressWarnings(sim <- get_p_values(fit, n_sim = 2))

  expect_identical(sim$simulation_fixef[[1]], c(NA, NA))

  expect_identical(
    sim$simulation_vcov[[1]],
    matrix(NA, nrow = 2, ncol = 2)
  )

  rm(
    list = c(
      "simulate.foo", "update.foo", "vcov.foo", "coef.foo", "formula.foo",
      "nobs.foo"
    ),
    envir = .GlobalEnv
  )
})
