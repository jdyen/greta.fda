#' A greta_fda regression model with function-valued data
#'
#' @description A \code{greta_fda} object contains a function regression model fitted with \code{greta_fda}
#' 
#' @rdname greta_fda
#' 
#' @param formula formula describing the model to be fitted, in the format used by \link[lme4]{lmer}
#' @param data a named list containing the variables in \code{formula}
#' @param family a GLM family passed as a \code{character}, see \link[stats]{family} (currently only gaussian (default) and poisson families are implemented)
#' @param y a \code{matrix} or \code{data.frame} with the response variable (one row per observation)
#' @param x a \code{matrix} or \code{data.frame} of predictor variables
#' @param z a \code{matrix} or \code{data.frame} of random effects variables
#' @param bins a vector of values at which the response variable is recorded
#' @param greta_settings a named list of values to pass to the inference method (see \link[greta]{inference} for details)
#' @param spline_settings a named list of settings to pass to the spline function (see \link[splines]{bs} for details)
#' @param priors a named list of prior distributions in the format of \link[greta]{distributions}
#' @param errors a character denoting the type of errors; currently only 'iid' is implemented
#' @param model a fitted \code{greta_fda} model
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{greta_fda}, which has associated `print`, `plot`, and `summary` methods
#' 
#' @export
#' 
#' @import greta
#' @importFrom splines bs
#' 
#' @examples
#' 
#' library(greta.fda)
#' 
#' # fit an example model
#' model <- greta_fda(y ~ x1 + x2 + (1 | z1),
#'                    data = greta_fda_data,
#'                    family = "poisson",
#'                    greta_settings = list(n_samples = 100,
#'                                          warmup = 100,
#'                                          chains = 2),
#'                    priors = list(intercept = normal(0, 1),
#'                                  coefficients = normal(0, 1),
#'                                  random_intercepts = normal(0, 1)))
#'                         
#' \dontrun{                 
#' # summarise fitted model
#' model
#' summary(model)
#' plot(model)
#' }

greta_fda <- function (y, ...) {
  
  UseMethod('greta_fda')
  
}

#' @rdname greta_fda
#'
#' @export
#' 
#' @examples
#'
#' # fit a greta_fda model with formula
#'   
#' \dontrun{
#' model <- greta_fda(y ~ x1 + x2 + (1 | z1),
#'                    data = greta_fda_data,
#'                    greta_settings = list(n_samples = 100,
#'                                          warmup = 100,
#'                                          chains = 2))
#' }

greta_fda.formula <- function (formula, data,
                               family = "gaussian",
                               bins = NULL,
                               greta_settings = list(),
                               spline_settings = list(),
                               priors = list(),
                               errors = 'iid',
                               ...) {
  
  # parse formula
  
  
  # create x, y, z objects to pass to default method
  y <- NULL
  x <- NULL
  z <- NULL

  # fit model
  model <- greta_fda.default(y = y, x = x, z = z,
                             family,
                             bins = bins,
                             greta_settings = greta_settings,
                             spline_settings = spline_settings,
                             priors = priors,
                             errors,
                             ...)
  
  # return fitted model
  model
  
  
}

#' @rdname greta_fda
#'
#' @export
#' 
#' @examples
#'
#' # fit a greta_fda model with separate arguments
#'   
#' \dontrun{
#' model <- greta_fda(y = y, x = x, z = z,
#'                    greta_settings = list(n_samples = 100,
#'                                          warmup = 100,
#'                                          chains = 2))
#' }

greta_fda.default <- function (y, x, z = NULL,
                               family = "gaussian",
                               bins = NULL,
                               greta_settings = list(),
                               spline_settings = list(),
                               priors = list(),
                               errors = 'iid', ...) {
  
  # test inputs (dimensions of y, x, z; classes of y, x, z)
  if (!all_equal(nrow(y), nrow(x), nrow(z))) {
    stop("x, y, and z must have the same number of rows")
  }
  if (!is.matrix(y)) {
    if (is.data.frame(y)) {
      classes <- apply(y, 2, class)
      if (!all(classes == "numeric")) {
        stop(paste0("the following columns in y are not numeric: ",
                    colnames(y)[which(classes != "numeric")]))
      }
      y <- as.matrix(y)
    } else {
      stop("y must be a matrix or data.frame")
    }
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x)) {
      classes <- apply(x, 2, class)
      if (!all(classes == "numeric")) {
        stop(paste0("x must be a valid design matrix but the following columns in x are not numeric: ",
                    colnames(x)[which(classes != "numeric")]))
      }
      x <- as.matrix(x)
    } else {
      if (length(x) == nrow(y)) {
        x <- matrix(x, ncol = 1)
      } else {
        stop("x must be a matrix or data.frame")
      }
    }
  }
  if (!is.null(z)) {
    if (!is.matrix(z)) {
      if (is.data.frame(z)) {
        classes <- apply(z, 2, class)
        if (!all(classes == "numeric")) {
          stop(paste0("z must be a valid design matrix but the following columns in z are not numeric: ",
                      colnames(z)[which(classes != "numeric")]))
        }
        z <- as.matrix(z)
      } else {
        if (length(z) == nrow(y)) {
          z <- matrix(z, ncol = 1)
        } else {
          stop("z must be a matrix or data.frame")
        }
      }
    }
    z <- apply(z, 2, function(x) as.integer(as.factor(z)))
  }

  # unpack greta_settings
  greta_set <- list(sampler = hmc(),
                    n_samples = 1000,
                    thin = 1,
                    warmup = 1000,
                    chains = 1,
                    verbose = TRUE,
                    pb_update = 50,
                    initial_values = NULL)
  greta_set[names(greta_settings)] <- greta_settings
  
  # unpack spline settings
  spline_set <- list(basis = bs,
                     df = 10,
                     degree = 3)
  spline_set[names(spline_settings)] <- spline_settings
  
  # prepare model
  greta_model <- build_greta_fda(y, x, z,
                                 family,
                                 bins,
                                 priors,
                                 errors,
                                 spline_set,
                                 ...)
  
  # sample from greta model
  samples <- mcmc(greta_model,
                  sampler = greta_set$sampler,
                  n_samples = greta_set$n_samples,
                  thin = greta_set$thin,
                  warmup = greta_set$warmup,
                  chains = greta_set$chains,
                  verbose = greta_set$verbose,
                  pb_update = greta_set$pb_update,
                  initial_values = greta_set$initial_values)
  
  # compile results
  model <- list(samples = samples,
                data = list(y = y,
                            x = x,
                            z = z),
                family = family,
                greta_settings = greta_set,
                spline_settings = spline_set,
                priors = priors,
                errors = errors)
  
  as.greta_fda(model)
  
}

#' @rdname greta_fda
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is a greta_fda model
#'   
#' \dontrun{
#' is.greta_fda(model)
#' }

is.greta_fda <- function (model) {
  inherits(x, 'greta_fda')
}

#' @rdname greta_fda
#'
#' @export
#'
#' @examples
#' 
#' # Print information about a 'greta_fda' object
#'
#' \dontrun{
#' print(model)
#' }

print.greta_fda <- function (model, ...) {
  cat(paste0('This is a greta_fda model\n'))
}

#' @rdname greta_fda
#'
#' @export
#'
#' @examples
#' 
#' # Plot a 'greta_fda' object
#'
#' \dontrun{
#' plot(model)
#' }

plot.greta_fda <- function (model, ...) {

  plot(model$samples, ...)

}


# internal function: create greta model from input data
build_greta_fda <- function (y, x, z,
                             family,
                             bins,
                             priors,
                             errors, 
                             spline_settings, ...) {
  
  # pull out index counters
  n <- nrow(y)
  nj <- ncol(y)
  nk <- ncol(x)
  if (!is.null(z)) {
    nt <- ncol(z)
    ngroup <- apply(z, 2, max)
  }

  # set up spline settings (nspline, nknots, degree)
  if (is.null(bins)) {
    bins <- seq_len(nj)
  }
  boundary_knots <- c(0, (nj + 1))
  np <- spline_settings$df

  # create spline basis
  spline_basis <- spline_settings$basis(bins,
                                        df = np,
                                        degree = splines_settings$degree,
                                        intercept = FALSE,
                                        Boundary.knots = boundary_knots)

  # setup priors
  sigma_main <- greta::uniform(min = 0.0, max = 5.0, dim = 1)
  sigma_bins <- greta::uniform(min = 0.0, max = 5.0, dim = nj)
  if (!is.null(z)) {
    sigma_gamma <- greta::uniform(min = 0.0, max = 5.0, dim = nt)
  }
  
  # setup parameters
  alpha <- greta::normal(mean = 0.0, sd = 1.0, dim = np)
  beta <- greta::normal(mean = 0.0, sd = 1.0, dim = c(nk, np))
  if (!is.null(z)) {
    gamma <- greta::normal(mean = rep(0.0, times = (nt * np)),
                           sd = rep(sigma_gamma, times = np),
                           dim = c(nt, np))
  }

  # define linear predictor
  mu <- greta::greta_array(0, dim = c(n, nj))
  for (i in seq_len(n)) {
    mu[i, ] <- sweep((x %*% (beta %*% spline_basis)), 2, (alpha %*% spline_basis), '+')
    if (!is.null(z)) {
      mu[i, ] <- mu[i, ] + (z %*% (gamma %*% spline_basis))
    }
  }
  
  # add error structure
  if (errors == 'iid') {
    bin_errors <- greta::normal(mean = rep(0.0, nj), sd = sigma_bins, dim = nj)
  } else {
    if (errors == 'ar1') {
      rho <- greta::uniform(min = 0, max = 1, dim = 1)
      bin_errors <- greta::normal(mean = rep(0.0, nj), sd = sigma_bins, dim = nj)
      bin_errors <- c(bin_errors[1], bin_errors[seq_len(nj)[-1]] + rho * bin_errors[seq_len(nj - 1)])
    } else {
      stop("errors must be 'iid' or 'ar1'")
    }
  }
  mu <- sweep(mu, 2, bin_errors, '+')
  
  # setup likelihood
  mu_vec <- do.call('c', mu)
  y_vec <- c(y)
  if (family == 'gaussian') {
    distribution(y_vec) <- greta::normal(mean = mu_vec, sd = sigma_main)
  } else {
    if (family == 'poisson') {
      distribution(y_vec) <- greta::poisson(mean = exp(mu_vec))
    } else {
      stop("family must be 'gaussian' or 'poisson'")
    }
  }
  
  # define model
  greta_model <- greta::model(alpha, beta, gamma, ...)

  # return model
  greta_model
  
}

# internal function: create greta_fda object
as.greta_fda <- function (model) {
  as_class(model, name = 'greta_fda', type = 'list')
}
