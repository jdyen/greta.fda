#' A fda_response object containing a greta_array for use with function-valued response variables
#'
#' @description A \code{fda_response} object contains a function regression object
#'   that can be used to model a function-valued response with \link[greta]{model}
#' 
#' @rdname fda_response
#' 
#' @param formula formula describing the model to be fitted, in the format used by \link[lme4]{lmer}
#' @param data a named list containing the variables in \code{formula}
#' @param y a \code{matrix} or \code{data.frame} with the response variable (one row per observation)
#' @param x a \code{matrix} or \code{data.frame} of predictor variables or a fda_response object for is.fda_response and print methods
#' @param z a \code{matrix} or \code{data.frame} of random effects variables
#' @param bins a vector of values at which the response variable is recorded
#' @param spline_settings a named list of settings to pass to the spline function (see \link[splines]{bs} for details)
#' @param priors a named list of settings for prior distributions (alpha_mean, alpha_sd, beta_mean, beta_sd, sigma_max)
#' @param errors a character denoting the type of errors in a matrix model; currently only 'iid' and 'ar1' are implemented
#' @param object an \code{fda_response} object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{fda_response}, which can be used in a \link[greta]{model}
#' 
#' @export
#' 
#' @import greta
#' @importFrom splines bs
#' @importFrom stats terms.formula delete.response as.formula model.matrix
#' 
#' @examples
#' 
#' library(greta.fda)
#' 
#' # fit an example model
#' fda_response <- fda_response(y ~ x1 + x2 + (1 | z1),
#'                              data = example_fda_data,
#'                              priors = list(alpha_mean = 0,
#'                                            alpha_sd = 1,
#'                                            beta_mean = 0,
#'                                            beta_sd = 1,
#'                                            sigma_max = 5))
#'                         
#' \dontrun{                 
#' # fit a greta model
#' sigma_response <- uniform(min = 0.0, max = 5.0, dim = 1)
#' distribution(example_fda_data$y) <- normal(fda_response$mu, sd = sigma_response)
#' 
#' greta_model <- with(fda_response, model(mu, alpha, beta, sigma_main, sigma_response))
#' }

fda_response <- function (y, ...) {
  
  UseMethod('fda_response')
  
}

#' @rdname fda_response
#'
#' @export
#' 
#' @examples
#'
#' # create a fda_response object to model a function-valued response variable
#'   
#' \dontrun{
#' 
#' # fit an example model
#' mu_fda <- fda_response(y ~ x1 + x2 + (1 | z1),
#'                        data = example_fda_data,
#'                        priors = list(alpha_mean = 0,
#'                                      alpha_sd = 1,
#'                                      beta_mean = 0,
#'                                      beta_sd = 1,
#'                                      sigma_max = 5))
#'                                      
#' sigma_response <- uniform(min = 0.0, max = 5.0, dim = 1)
#' distribution(example_fda_data$y) <- normal(mu_fda, sd = sigma_response)
#' 
#' greta_model <- model(mu_fda)
#' 
#' }

fda_response.formula <- function (formula, data,
                                  bins = NULL,
                                  spline_settings = list(),
                                  priors = list(),
                                  errors = 'iid',
                                  ...) {
  
  # parse formula
  response <- all.vars(formula)[1] 
  terms <- terms(formula)
  random <- (grep("\\|", attributes(terms)$term.labels))
  var_names <- all.vars(delete.response(terms))
  
  # use correct var_names when random is missing
  if (length(random)) {
    fixed_vars <- var_names[-random]
    random_vars <- var_names[random]
  } else {
    fixed_vars <- var_names
    random_vars <- NULL
  }

  # create x, y, z objects to pass to default method
  y <- get(response, envir = as.environment(data), inherits = TRUE)
  if (length(fixed_vars)) {
    x_tmp <- mget(fixed_vars, envir = as.environment(data), inherits = TRUE)
  }
  if (length(random_vars)) {
    z_tmp <- mget(random_vars, envir = as.environment(data), inherits = TRUE)
    z_tmp <- lapply(z_tmp, function(x) as.integer(as.factor(x)))
  }
  
  # create model matrix
  if (length(fixed_vars)) {
    x <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(fixed_vars, collapse = " + "))), data = x_tmp)
  } else {
    x <- matrix(0, nrow = length(y), ncol = 1)
    colnames(x) <- 'null'
  }
  if (length(random_vars)) {
    z <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(random_vars, collapse = " + "))), data = z_tmp)
  } else {
    z <- NULL
  }
  
  # create fda_response object
  fda_response <- fda_response.default(y = y, x = x, z = z,
                                       bins = bins,
                                       spline_settings = spline_settings,
                                       priors = priors,
                                       errors = errors,
                                       ...)
  
  # add formula to fda_response object
  fda_response$formula <- formula
  
  # return fda_response object
  fda_response
  
}

#' @rdname fda_response
#'
#' @export
#' 
#' @examples
#'
#' # default method to create a fda_response object
#'   
#' \dontrun{
#' 
#' # fit an example model
#' y <- example_fda_data$y
#' xmat <- cbind(example_fda_data$x1, example_fda_data$x2)
#' xmat <- example_fda_data$z1
#' mu_fda <- fda_response(y = y, x = xmat, z = zmat,
#'                        priors = list(alpha_mean = 0,
#'                                      alpha_sd = 1,
#'                                      beta_mean = 0,
#'                                      beta_sd = 1,
#'                                      sigma_max = 5))
#'                                      
#' sigma_response <- uniform(min = 0.0, max = 5.0, dim = 1)
#' distribution(example_fda_data$y) <- normal(mu_fda, sd = sigma_response)
#' 
#' greta_model <- model(mu_fda)
#' 
#' }

fda_response.default <- function (y, x, z = NULL,
                                  bins = NULL,
                                  spline_settings = list(),
                                  priors = list(),
                                  errors = 'iid', ...) {
  
  # test inputs (dimensions of y, x, z; classes of y, x, z)
  if (is.null(z)) {
    nrow_z <- nrow(x)
  } else { 
    nrow_z <- nrow(z)
  } 
  if (!is.matrix(y)) {
    if (is.data.frame(y)) {
      classes <- apply(y, 2, class)
      if (!all(classes == "numeric")) {
        stop(paste0('the following columns in y are not numeric: ',
                    colnames(y)[which(classes != 'numeric')]))
      }
      y <- as.matrix(y)
      if (!all_equal(nrow(y), nrow(x), nrow_z)) {
        stop('x, y, and z must have the same number of rows')
      }
      model_type <- 'matrix'
    } else {
      if (is.numeric(y)) {
        if (all_equal(length(y), nrow(x), nrow_z)) {
          model_type <- 'flat'
          if (is.null(bins)) {
            stop('bins must be provided if the response variable is flattened')
          }
        } else {
          stop('length of y must be equal to the number of rows in x and z')
        }
      } else {
        stop('y must be a numeric matrix, vector, or data.frame')
      }
    }
  } else {
    if (!all_equal(nrow(y), nrow(x), nrow_z)) {
      stop("x, y, and z must have the same number of rows")
    }
    model_type <- 'matrix'
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
    z <- apply(z, 2, function(x) as.integer(as.factor(x)))
  }

  # unpack spline settings
  spline_set <- list(basis = "bs",
                     df = 10,
                     degree = 3)
  spline_set[names(spline_settings)] <- spline_settings
  
  # prepare greta_array
  if (model_type == 'matrix') {
    fda_response <- build_fda_response_matrix(y, x, z,
                                              bins,
                                              priors,
                                              errors,
                                              spline_set,
                                              ...)
  }
  if (model_type == 'flat') {
    fda_response <- build_fda_response_flat(y, x, z,
                                            bins,
                                            priors,
                                            errors,
                                            spline_set,
                                            ...)
  }  
  
  # add variable names to fda_response object
  fda_response$var_names <- list(x = colnames(x),
                                 z = colnames(z))
  
  as.fda_response(fda_response)
  
}

#' @rdname fda_response
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is a fda_response object
#'   
#' \dontrun{
#' is.fda_response(x)
#' }

is.fda_response <- function (x) {
  inherits(x, 'fda_response')
}

#' @rdname fda_response
#'
#' @export
#'
#' @examples
#' 
#' # Print information about a 'fda_response' object
#'
#' \dontrun{
#' print(x)
#' }

print.fda_response <- function (x, ...) {
  cat(paste0('This is a fda_response object\n'))
}

#' @rdname fda_response
#'
#' @export
#'
#' @examples
#' 
#' # Summarise a 'fda_response' object
#'
#' \dontrun{
#' summary(object)
#' }

summary.fda_response <- function (object, ...) {
  cat(paste0('This is a fda_response object\n'))
}


# internal function: create greta array from matrix input data
build_fda_response_matrix <- function (y, x, z,
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
  
  # convert x and y to greta arrays
  x <- greta::as_data(x)
  y <- greta::as_data(y)
  
  # set up spline settings (nspline, nknots, degree)
  if (is.null(bins)) {
    bins <- seq_len(nj)
    max_bound <- nj + 1
  } else {
    max_bound <- max(bins) + 1
  }
  boundary_knots <- c(0, max_bound)
  np <- spline_settings$df

  # create spline basis
  spline_basis <- get(spline_settings$basis)(bins,
                                             df = np,
                                             degree = spline_settings$degree,
                                             intercept = FALSE,
                                             Boundary.knots = boundary_knots)
  spline_basis <- greta::as_data(t(spline_basis))

  # setup priors
  prior_set <- list(alpha_mean = 0.0,
                    alpha_sd = 1.0,
                    beta_mean = 0.0,
                    beta_sd = 1.0,
                    sigma_max = 5.0)
  prior_set[names(priors)] <- priors
  sigma_main <- greta::uniform(min = 0.0, max = prior_set$sigma_max, dim = 1)
  sigma_bins <- greta::uniform(min = 0.0, max = prior_set$sigma_max, dim = nj)
  sigma_gamma <- NULL
  if (!is.null(z)) {
    sigma_gamma <- greta::uniform(min = 0.0, max = prior_set$sigma_max, dim = c(nt, np))
  }
  
  # setup parameters
  alpha <- greta::normal(mean = prior_set$alpha_mean, sd = prior_set$alpha_sd, dim = c(1, np))
  beta <- greta::normal(mean = prior_set$beta_mean, sd = prior_set$beta_sd, dim = c(nk, np))
  
  if (!is.null(z)) {
    gamma <- vector('list', length = nt)
    for (rand in seq_len(nt)) {
      gamma[[rand]] <- greta::normal(mean = greta::zeros(dim = c(ngroup[rand], np)),
                                     sd = greta::greta_array(rep(sigma_gamma[rand, ], ngroup[rand]),
                                                             dim = c(ngroup[rand], np)),
                                     dim = c(ngroup[rand], np))
    }
  }
  
  # define linear predictor
  mu <- sweep((x %*% (beta %*% spline_basis)), 2, t(alpha %*% spline_basis), '+')
  if (!is.null(z)) {
    for (rand in seq_len(nt)) {
      mu <- mu + (gamma[[rand]][z[, rand], ] %*% spline_basis)
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
  
  # flatten gamma list (if used)
  gamma_vec <- NULL
  if (!is.null(z)) {
    gamma_vec <- do.call('c', gamma)
  }

  # return model
  list(mu = mu,
       alpha = alpha,
       beta = beta,
       gamma = gamma_vec,
       sigma_main = sigma_main,
       sigma_bins = sigma_bins,
       sigma_gamma = sigma_gamma,
       bins = bins,
       spline_basis = spline_basis,
       spline_settings = spline_settings,
       data = list(y = y,
                   x = x,
                   z = z),
       formula = NULL,
       priors = prior_set, 
       errors = errors)
  
}

# internal function: create greta array from flattened input data
build_fda_response_flat <- function (y, x, z,
                                     bins,
                                     priors,
                                     errors, 
                                     spline_settings, ...) {
  
  # pull out index counters
  n <- length(y)
  nk <- ncol(x)
  if (!is.null(z)) {
    nt <- ncol(z)
    ngroup <- apply(z, 2, max)
  }
  
  # convert x and y to greta arrays
  x <- greta::as_data(x)
  y <- greta::as_data(y)

  # set up spline settings (nspline, nknots, degree)
  boundary_knots <- c(0, max(bins) + 1)
  np <- spline_settings$df
  
  # create spline basis
  spline_basis <- get(spline_settings$basis)(bins,
                                             df = np,
                                             degree = spline_settings$degree,
                                             intercept = FALSE,
                                             Boundary.knots = boundary_knots)
  spline_basis <- greta::as_data(t(spline_basis))

  # setup priors
  prior_set <- list(alpha_mean = 0.0,
                    alpha_sd = 1.0,
                    beta_mean = 0.0,
                    beta_sd = 1.0,
                    sigma_max = 5.0)
  prior_set[names(priors)] <- priors
  sigma_main <- greta::uniform(min = 0.0, max = prior_set$sigma_max, dim = 1)
  sigma_bins <- NULL
  sigma_gamma <- NULL
  if (!is.null(z)) {
    sigma_gamma <- greta::uniform(min = 0.0, max = prior_set$sigma_max, dim = c(nt, np))
  }
  
  # setup parameters
  alpha <- greta::normal(mean = prior_set$alpha_mean, sd = prior_set$alpha_sd, dim = c(1, np))
  beta <- greta::normal(mean = prior_set$beta_mean, sd = prior_set$beta_sd, dim = c(nk, np))
  
  if (!is.null(z)) {
    gamma <- vector('list', length = nt)
    for (rand in seq_len(nt)) {
      gamma[[rand]] <- greta::normal(mean = greta::zeros(dim = c(ngroup[rand], np)),
                                     sd = greta::greta_array(rep(sigma_gamma[rand, ], ngroup[rand]),
                                                             dim = c(ngroup[rand], np)),
                                     dim = c(ngroup[rand], np))
    }
  }

  # define linear predictor
  mu <- t(alpha %*% spline_basis) + rowSums(x * t(beta %*% spline_basis))
  if (!is.null(z)) {
    for (rand in seq_len(nt)) {
      mu <- mu + rowSums(gamma[[rand]][z[, rand], ] * t(spline_basis))
    }
  }
  
  # flatten gamma list (if used)
  gamma_vec <- NULL
  if (!is.null(z)) {
    gamma_vec <- do.call('c', gamma)
  }
  
  # return model
  list(mu = mu,
       alpha = alpha,
       beta = beta,
       gamma = gamma_vec,
       sigma_main = sigma_main,
       sigma_bins = sigma_bins,
       sigma_gamma = sigma_gamma,
       bins = bins,
       spline_basis = spline_basis,
       spline_settings = spline_settings,
       data = list(y = y,
                   x = x,
                   z = z),
       formula = NULL,
       priors = prior_set, 
       errors = errors)
  
}

# internal function: create a fda_response object
as.fda_response <- function (x) {
  as_class(x, name = 'fda_response', type = 'list')
}
