% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_refit.R
\name{get_refit}
\alias{get_refit}
\alias{get_refit.default}
\alias{get_refit.merMod}
\alias{get_refit.glmmTMB}
\title{Refit a model}
\usage{
get_refit(object, new_response, ...)

\method{get_refit}{default}(object, new_response, ...)

\method{get_refit}{merMod}(object, new_response, ...)

\method{get_refit}{glmmTMB}(object, new_response, ...)
}
\arguments{
\item{object}{A model.}

\item{new_response}{the new response, may be a vector or a matrix.}

\item{...}{other arguments passed to \code{refit} or \code{update}.}
}
\value{
A model with same class as \code{object}.
}
\description{
Refit a model with a new response.
}
\details{
This function uses \code{new_response} to refit \code{object} replacing its old response variable.
If the class is \code{merMod} it uses \code{refit}, otherwise uses \code{\link[stats:update]{stats::update()}}.

The default method tries to update the model response using it's \code{\link[stats:model.frame]{stats::model.frame()}},
if it errors it tries to update the model by inserting the \code{new_response}
directly into the object formula.
}
\seealso{
\code{\link[stats:update]{stats::update()}}
}
