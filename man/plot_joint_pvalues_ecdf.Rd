% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_size.R
\name{plot_joint_pvalues_ecdf}
\alias{plot_joint_pvalues_ecdf}
\title{Plot Simulation joint p-values}
\usage{
plot_joint_pvalues_ecdf(
  model,
  generator = NULL,
  n_sim = 1000,
  ylab = "Fn(x)",
  xlab = "x",
  ...
)
}
\arguments{
\item{model}{A model with [update()] method.}

\item{generator}{Optional custom generator function with 2 arguments (n, mu).}

\item{n_sim}{Number of simulations.}

\item{ylab}{y-axis label}

\item{xlab}{x-axis label}

\item{...}{Extra arguments from [plot()]}
}
\value{
Double vector with joint p-values.
}
\description{
Plot Simulation joint p-values
}
\examples{
fit <- lm(mpg ~ cyl, data = mtcars)

plot_joint_pvalues_ecdf(fit)
}
