% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fda_predictor-class.R
\name{fda_predictor}
\alias{fda_predictor}
\alias{is.fda_predictor}
\alias{print.fda_predictor}
\alias{summary.fda_predictor}
\title{A fda_predictor object containing a greta_array for use with function-valued predictor variables}
\usage{
fda_predictor(x, bins = NULL, spline_settings = list(), priors = list(), ...)

is.fda_predictor(x)

\method{print}{fda_predictor}(x, ...)

\method{summary}{fda_predictor}(object, ...)
}
\arguments{
\item{x}{a \code{numeric} vector, \code{matrix}, or \code{data.frame} with a single predictor variable (or a fda_predictor object for is.fda_predictor and print methods)}

\item{bins}{a \code{integer} or \code{numeric} vector or \code{matrix} of values at which the predictor variable is recorded}

\item{spline_settings}{a named list of settings to pass to the spline function (see \link[splines]{bs} for details)}

\item{priors}{a named list of settings for prior distributions (mean, sd, sigma_max)}

\item{...}{further arguments passed to or from other methods}

\item{object}{a \code{fda_predictor} object}
}
\value{
An object of class \code{fda_predictor}, which can be used in a \link[greta]{model}
}
\description{
A \code{fda_predictor} object contains a function regression object
  that can be used to model a function-valued predictor with \link[greta]{model}
}
\examples{

library(greta.fda)

# fit an example model
fda_x <- fda_predictor(example_fda_data$x_fun,
                             priors = list(
                               nu_local = 2,
                               nu_global = 2,
                               slab_df = 2,
                               scale_global = 2,
                               slab_scale = 2
                             )))
                        
\dontrun{                 

# define some priors
alpha <- normal(mean = 0.0, sd = 1.0, dim = 1)
beta <- normal(mean = 0.0, sd = 1.0, dim = 1)
sigma_response <- uniform(min = 0.0, max = 5.0, dim = 1)

# define a linear predictor with the fda_predictor object
mu <- alpha + beta * example_fda_data$x1 + fda_x$array

# fit a model
distribution(example_fda_data$y_scalar) <- normal(mu, sd = sigma_response)

greta_model <- model(mu, alpha, beta, fda_x$beta, sigma_response)
}

# check if an object is a fda_predictor object
  
\dontrun{
is.fda_predictor(x)
}

# Print information about a 'fda_predictor' object

\dontrun{
print(x)
}

# Summarise a 'fda_predictor' object

\dontrun{
summary(object)
}
}
