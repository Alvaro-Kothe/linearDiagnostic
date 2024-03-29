% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_size.R
\name{get_p_values_joint}
\alias{get_p_values_joint}
\title{Generate P-Values from Model Coefficient Simulations}
\usage{
get_p_values_joint(model, n_sim = 1000, test_coefficients = NULL, ...)
}
\arguments{
\item{model}{A model with \code{\link[=update]{update()}}, \code{\link[=simulate]{simulate()}}, \code{\link[=fitted]{fitted()}},
\code{\link[=coef]{coef()}}, \code{\link[=vcov]{vcov()}} methods.}

\item{n_sim}{The number of simulations to perform.}

\item{test_coefficients}{Numeric vector. A vector with values to be used to compute
the test statistic. It should be the coefficients that was used to compute
the fitted values of the response. If \code{NULL} defaults to coef(model)}

\item{...}{Additional arguments to be passed to \code{simulate_coefficients}.}
}
\value{
A numeric vector containing the joint p-values obtained from each simulation.
}
\description{
This function generates p-values by simulating coefficients from
a given model and computing Wald statistics for each simulation.
}
\examples{
\dontrun{
model <- lm(mpg ~ wt + hp, data = mtcars)
get_p_values_joint(model, n_sim = 100)
}

}
