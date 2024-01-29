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
