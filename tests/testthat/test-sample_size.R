test_that("simulate_coefficients() returns a list of size `n_sim`", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)

  simulation_result <- simulate_coefficients(fit, n_sim = 5)

  expect_named(simulation_result, c("coefs", "vcov"))
  expect_length(simulation_result$coefs, 5)
  expect_length(simulation_result$vcov, 5)
})

test_that("simulate_coefficients() coefs length is equal to number of parameters estimated", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)

  simulation_result <- simulate_coefficients(fit, n_sim = 1)

  expect_length(simulation_result$coefs[[1]], 2)
  expect_identical(dim(simulation_result$vcov[[1]]), c(2L, 2L))
})

test_that("set.seed works with simulate_coefficients()", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)

  set.seed(1)
  simulation_result1 <- simulate_coefficients(fit, n_sim = 2)

  set.seed(1)
  simulation_result2 <- simulate_coefficients(fit, n_sim = 2)

  expect_identical(simulation_result1, simulation_result2)
})

test_that("simulate_coefficients() with generator yield different results", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)

  set.seed(1)
  simulation_result1 <- simulate_coefficients(fit, n_sim = 2)

  set.seed(1)
  simulation_result2 <- simulate_coefficients(fit,
    n_sim = 2,
    generator = function(n, mu) rgamma(n, mu)
  )

  expect_false(identical(simulation_result1, simulation_result2))
})

test_that("plot_*_pvalues_ecdf work", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)
  mpg_fit <- lm(mpg ~ ., data = mtcars)

  set.seed(1)
  expect_no_error(plot_joint_pvalues_ecdf(fit, n_sim = 10))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10, ask = TRUE))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10, ask = FALSE))
  expect_no_error(plot_joint_pvalues_ecdf(fit, n_sim = 10, plot_uniform = TRUE, uniform_legend = FALSE))
  expect_no_error(plot_joint_pvalues_ecdf(fit, n_sim = 10, plot_uniform = TRUE))
  expect_no_error(plot_joint_pvalues_ecdf(fit, n_sim = 10, plot_uniform = FALSE))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10, plot_uniform = TRUE))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10, plot_uniform = TRUE, uniform_legend = FALSE))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10, plot_uniform = FALSE))
})

test_that("plot_*_pvalues throw warning with singular matrix", {
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
      plot_joint_pvalues_ecdf(fit, n_sim = 10),
      "Couldn't inverse vcov from \\d+ simulations and used ginv instead"
    )
  })
})

test_that("simulate_coefficients() works with poisson with offset", {
  n <- 5
  offset_ <- rpois(n, 10000)
  x <- runif(n)
  y <- rpois(n, exp(1 + 2 * x))

  suppressWarnings({
    fit <- glm(y ~ x + offset(offset_), family = poisson())
    sim <- simulate_coefficients(fit, n_sim = 1)
  })

  expect_length(sim$coefs[[1]], 2)
})

test_that("plot_*_ecdf() errors when data is unavailable", {
  df <- data.frame(
    x = c(1, 3, 5, 9),
    y = c(1, 2, 3, 4)
  )
  fit <- lm(y ~ x, data = df)
  rm(df)

  expect_error(plot_joint_pvalues_ecdf(fit, n_sim = 10))
  expect_error(plot_pvalues_ecdf(fit, n_sim = 10))
  expect_no_error(plot_joint_pvalues_ecdf(fit, n_sim = 10, args_sim = list(data = model.frame(fit))))
  expect_no_error(plot_pvalues_ecdf(fit, n_sim = 10, args_sim = list(data = model.frame(fit))))
})


test_that("simulate_coefficients() with custom method works", {
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

    class(out) <- "foo"

    return(out)
  }

  # nolint start: object_name_linter
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
  # nolint end


  fit <- foo(mpg ~ cyl, mtcars)

  sim <- simulate_coefficients(fit, n_sim = 2)

  expect_identical(sim$coefs, list(c(NA, NA), c(NA, NA)))

  expect_identical(sim$vcov, list(
    matrix(NA, nrow = 2, ncol = 2),
    matrix(NA, nrow = 2, ncol = 2)
  ))

  rm(
    list = c("simulate.foo", "update.foo", "vcov.foo", "coef.foo"),
    envir = .GlobalEnv
  )
})
