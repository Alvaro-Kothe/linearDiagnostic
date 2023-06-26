plot_cook <- function(model, n_highlights = 0, cut = FALSE,
                      xlab = "Index", ylab = "Cook's distance") {
  cook <- cooks.distance(model)

  plot(seq_along(cook), cook, pch = 20, xlab = xlab, ylab = ylab)
  if (n_highlights > 0) {
    highests <- order(cook, decreasing = TRUE)[seq_len(n_highlights)]
    text(x = highests, y = cook[highests], label = highests)
  }
  if (cut) {
    abline(h = mean(cook) + 4 * sd(cook))
  }
}

plot_res_vs_linear_predictor <- function(model, residual_fn = rstandard,
                                         type = "deviance",
                                         xlab = "Linear Predictor",
                                         ylab = "Standardized deviance residuals",
                                         main = "") {
  linear_predictor <- switch(class(model),
    glm = predict(model, type = "link"),
    fitted(model)
  )

  plot(x = linear_predictor, residual_fn(model))
}
