test_that("Plot functions runs without problems", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)

  expect_no_error(plot_cook(fit))
  expect_no_error(plot_cook(fit, n_highlights = length(x) + 1))
  expect_no_error(plot_cook(fit, n_highlights = 2))
  expect_no_error(plot_cook(fit, cut = TRUE))
  expect_no_error(plot_res_vs_linear_predictor(fit))
  expect_no_error(plot_res_vs_linear_predictor(fit, pch = 20))
})

test_that("plot_*_pvalues_ecdf work", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 3, 6, 9)
  fit <- lm(y ~ x)
  p_values <- get_p_values(fit, n_sim = 10)
  expect_no_error(plot(p_values))
})
