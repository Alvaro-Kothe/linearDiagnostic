# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
envelope_measures <- function(fit.model, residual = "studentized",
                                type = "pearson", standardized = TRUE,
                                alpha = .05) {
  resid_fn <- switch(residual,
                     "standardized" = rstandard,
                     "studentized" = rstudent,
                     "regular" = resid
  )
  td <- resid_fn(fit.model, type = type)
  model_frame <- fit.model$model
  n <- length(td)
  e <- matrix(0, n, 100)
  #
  for (i in 1:100) {
    y_ <- simulate(fit.model)
    model_frame[, 1] <- y_
    fit <- update(fit.model, formula = . ~ ., data = model_frame)
    e[, i] <- sort(resid_fn(fit, type = type))
  }
  #
  es <- apply(e, 1, function(x) quantile(x, c(alpha / 2, .5, 1 - alpha / 2)))

  tdsort <- sort(td)
  oob <- !(tdsort > es[1, ] & tdsort < es[3, ])

  qq <- qqnorm(td, plot.it = FALSE)
  qq_ord <- order(qq$x)
  qq_sort <- qq$x[qq_ord]

  return(list(
    expected = qq_sort, observed = qq$y[qq_ord], outside = oob,
    lower = es[1, ], med = es[2, ], upper = es[3, ]
  ))
}

# plot_envelope_ <- function(expected, observed, li, med, ls, color = NULL,
#                              xlab = "Expected quantiles",
#                              ylab = "Observed quantiles") {
#   out <- ggplot2::ggplot() +
#     ggplot2::geom_point(aes(expected, observed,
#                             colour = color
#     )) +
#     ggplot2::geom_line(aes(expected, li), linetype = 1) +
#     ggplot2::geom_line(aes(expected, med), linetype = 2) +
#     ggplot2::geom_line(aes(expected, ls), linetype = 1)
#
#   return(out)
# }

plot_envelope_base <- function(expected, observed, lower, med, upper,
                             xlab = "Expected quantiles",
                             ylab = "Observed quantiles") {

  plot(expected, observed, type = "p", pch = 20, xlab=xlab, ylab=ylab)
  lines(expected, lower, lty = 1)
  lines(expected, upper, lty = 1)
  lines(expected, med, lty = 2)
  return(invisible(list(x = expected, y = observed, lower=lower, upper=upper, med=med)))
}

plot_envelope <- function(model, xlab = "Expected quantiles", ylab = "Observed quantiles",
                          seed = NULL, ...) {

  set.seed(seed)
  env_meas <- envelope_measures(model, ...)

  envelope_plot <- plot_envelope_base(env_meas$expected, env_meas$observed,
                                 env_meas$lower, env_meas$med, env_meas$upper,
                                 xlab = xlab, ylab = ylab)


  return(envelope_plot)
}
