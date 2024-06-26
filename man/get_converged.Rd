% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_coefs_and_vcov.R
\name{get_converged}
\alias{get_converged}
\alias{get_converged.default}
\alias{get_converged.merMod}
\title{Did the model converged}
\usage{
get_converged(object)

\method{get_converged}{default}(object)

\method{get_converged}{merMod}(object)
}
\arguments{
\item{object}{A model to check for convergence.}
}
\value{
A boolean value.
}
\description{
Retrieves if the model converged.
}
\details{
The default method retrieves \code{object$converged}.
For models of class merMod it verifies if the infinity norm of the
Newthon-Raphson step is less than 0.0001.
}
