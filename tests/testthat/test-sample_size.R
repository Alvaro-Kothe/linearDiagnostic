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
