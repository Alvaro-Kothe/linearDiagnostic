% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/envelope.R
\name{envelope}
\alias{envelope}
\title{Generate Simulated Envelope}
\usage{
envelope(
  model,
  residual_fn = stats::rstudent,
  alpha = 0.05,
  n_sim = 100,
  plot.it = TRUE,
  ...
)
}
\arguments{
\item{model}{A model that has a \code{\link[stats:simulate]{stats::simulate()}} method and is compatible with \code{\link[=get_refit]{get_refit()}}.}

\item{residual_fn}{A function to calculate model residuals. The default is
\code{\link[stats:influence.measures]{stats::rstudent()}} for studentized residuals.}

\item{alpha}{The significance level for constructing the envelope bounds.
Defaults to 0.05.}

\item{n_sim}{The number of simulations to perform for envelope construction.
Defaults to 100.}

\item{plot.it}{Logical. Generate envelope plot.}

\item{...}{Extra arguments to \code{\link[=get_refit]{get_refit()}}}
}
\value{
An object of class \code{LD_envelope}, which contains the following components:
\describe{
\item{expected}{A vector of expected quantiles from a normal distribution.}
\item{observed}{A vector of observed quantiles from the model residuals.}
\item{outside}{A logical vector indicating whether each observation falls
outside the constructed envelope bounds.}
\item{lower}{The lower bounds of the envelope for each observation.}
\item{med}{The median bounds of the envelope for each observation.}
\item{upper}{The upper bounds of the envelope for each observation.}
}
}
\description{
Generates a normal QQ-plot with simulated envelope residuals.
}
\details{
Simulates new responses using \code{\link[stats:simulate]{stats::simulate()}} and refits the model
for each vector of new responses using \code{\link[=get_refit]{get_refit()}}. The function then computes
residuals for each simulation, sorts them, and constructs envelope bands and
a median line based on the quantiles of these residuals.
}
\examples{
\dontrun{
fit <- lm(mpg ~ cyl, data = mtcars)

envelope(fit)
envelope(fit, residual_fn = function(x) residuals.lm(x, type = "deviance"))
envelope(fit, residual_fn = rstandard, n_sim = 200)
}

}
\seealso{
\code{\link{get_refit}}, \code{\link{simulate}}, \code{\link{rstudent}}, \code{\link{plot.LD_envelope}}
}
