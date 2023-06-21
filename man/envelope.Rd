\name{plot_envelope}
\alias{plot_envelope}}
\title{Plot Envelope}}}
\usage{
plot_envelope(model, xlab = "Expected quantiles", ylab = "Observed quantiles",
                          seed = NULL, ...)
}
\description{
Prints 'Hello, world!'.
}
\examples{
fit <- lm(mpg~cyl, data = mtcars)
plot_envelope(fit)
}
