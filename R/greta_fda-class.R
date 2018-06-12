#' A greta_fda regression model with function-valued data
#'
#' @description A \code{greta_fda} object contains a function regression model fitted with \code{greta_fda}
#' 
#' @rdname greta_fda
#' 
#' @param formula formula describing the model to be fitted, in the format used by \link[lme4]{lmer}
#' @param data a named list containing the variables in \code{formula}
#' @param family a GLM family passed as a \code{character}, see \link[stats]{family} (currently only gaussian (default) and poisson families are implemented)
#' @param link a link function for the specified GLM family, see \link[stats]{family} for details
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
#' @importFrom stats terms.formula delete.response
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
                               family = 'gaussian',
                               link = NULL,
                               bins = NULL,
                               greta_settings = list(),
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
  
  # fit model
  model <- greta_fda.default(y = y, x = x, z = z,
                             family = family,
                             link = link,
                             bins = bins,
                             greta_settings = greta_settings,
                             spline_settings = spline_settings,
                             priors = priors,
                             errors,
                             ...)
  
  # add formula to fitted model
  model$formula <- formula
  
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
                               link = NULL,
                               bins = NULL,
                               greta_settings = list(),
                               spline_settings = list(),
                               priors = list(),
                               errors = 'iid', ...) {
  
  # test inputs (dimensions of y, x, z; classes of y, x, z)
  if (is.null(z)) {
    nrow_z <- nrow(x)
  } else {
    nrow_z <- nrow(z)
  }
  if (!all_equal(nrow(y), nrow(x), nrow_z)) {
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
    z <- apply(z, 2, function(x) as.integer(as.factor(x)))
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
  spline_set <- list(basis = "bs",
                     df = 10,
                     degree = 3)
  spline_set[names(spline_settings)] <- spline_settings
  
  # prepare model
  greta_model <- build_greta_fda(y, x, z,
                                 family,
                                 link,
                                 bins,
                                 priors,
                                 errors,
                                 spline_set,
                                 ...)
  
  # set link function (to be added)
  if (is.null(link)) {
    if (family == 'gaussian') {
      link <- 'identity'
    }
    if (family == 'poisson') {
      link <- 'log'
    }
    if (family == 'binomial') {
      link <- 'cloglog'
    }
  }
  
  # sample from greta model
  samples <- mcmc(greta_model$greta_model,
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
                            z = z,
                            bins = greta_model$bins),
                family = family,
                spline_basis = greta_model$spline_basis,
                greta_settings = greta_set,
                spline_settings = spline_set,
                link = link,
                formula = NULL,
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
  inherits(model, 'greta_fda')
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
#' print(x)
#' }

print.greta_fda <- function (x, ...) {
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
#' plot(x)
#' }

plot.greta_fda <- function (x, ...) {

  plot(x$samples, ...)

}

#' @rdname greta_fda
#'
#' @export
#'
#' @examples
#' 
#' # Summarise a 'greta_fda' object
#'
#' \dontrun{
#' summary(x)
#' }

summary.greta_fda <- function (x, ...) {
  
  # calculate Bayesian R2
  # y <- x$data$y
  # yhat <- predict(x)
  # e <- -1 * sweep(yhat, 2, y)
  # var_yhat <- apply(yhat, 1, var)
  # var_e <- apply(e, 1, var)
  # bayes_R2 <- var_yhat / (var_yhat + var_e)

  # calculate rhat
  rhat <- NULL
  
  # return outputs
  list(bayes_R2 = bayes_R2,
       coefficients = coef(x),
       rhat = rhat)
  
}

#' @rdname greta_fda
#'
#' @export
#'
#' @examples
#' 
#' # Extract coefficients from a 'greta_fda' object
#'
#' \dontrun{
#' coef(x)
#' }

coef.greta_fda <- function (x, ...) {
  
  # extract coefficients
  param_estimates <- lapply(x$samples, function(z) apply(z, 2, mean))
  beta <- lapply(param_estimates, function(z) matrix(z[grep('beta', names(z))],
                                                     nrow = ncol(x$data$x)))
  beta_mean <- array(NA, dim = c(nrow(beta[[1]]), ncol(beta[[1]]), length(beta)))
  for (i in seq_along(beta)) {
    beta_mean[, , i] <- beta[[i]]
  }
  beta_mean <- apply(beta_mean, c(1, 2), mean)
  beta_link <- beta_mean %*% as.matrix(x$spline_basis)
  
  # return some summary of this
  list(beta = beta_link)
  
}

#' @rdname greta_fda
#'
#' @export
#'
#' @examples
#' 
#' # Calculate fitted values from a 'greta_fda' object
#'
#' \dontrun{
#' fitted(x)
#' }

fitted.greta_fda <- function (x, ...) {
  
  # calculate fitted values
  param_estimates <- lapply(x$samples, function(z) apply(z, 2, mean))
  fitted <- lapply(param_estimates,
                   function(z) z[grep('mu', names(param_estimates[[1]]))])
  fitted <- do.call('rbind', fitted)
  fitted <- apply(fitted, 2, mean)
  
  # return some summary of fitted values
  matrix(fitted, ncol = ncol(x$data$y))
  
}

#' @rdname greta_fda
#'
#' @export
#'
#' @examples
#' 
#' # Predict values from a 'greta_fda' object
#'
#' \dontrun{
#' predict(x)
#' }

predict.greta_fda <- function (x, ..., newdata = NULL, type = c("link", "response"),
                               re.form = NULL, fun = NULL) {
  
  # fill data (check lme4 method for this)
  if (is.null(newdata)) {
    newdata <- x$data
  }
  
  # calculate coefficient
  coefs <- coef(x)
  
  # predict outputs
  out <- newdata$x %*% coefs$beta
  
  # add random effects based on re.form
  
  # set link function
  
  
  # return outputs
  out
  
}


# internal function: create greta model from input data
build_greta_fda <- function (y, x, z,
                             family,
                             link,
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
  spline_basis <- get(spline_settings$basis)(bins,
                                             df = np,
                                             degree = spline_settings$degree,
                                             intercept = FALSE,
                                             Boundary.knots = boundary_knots)
  spline_basis <- t(spline_basis)
  
  # setup priors
  sigma_main <- greta::uniform(min = 0.0, max = 5.0, dim = 1)
  sigma_bins <- greta::uniform(min = 0.0, max = 5.0, dim = nj)
  if (!is.null(z)) {
    sigma_gamma <- greta::uniform(min = 0.0, max = 5.0, dim = c(nt, np))
  }
  
  # setup parameters
  alpha <- greta::normal(mean = 0.0, sd = 1.0, dim = c(1, np))
  beta <- greta::normal(mean = 0.0, sd = 1.0, dim = c(nk, np))
  
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
  
  # setup likelihood
  if (family == 'gaussian') {
    distribution(y) <- greta::normal(mean = mu, sd = sigma_main)
  } else {
    if (family == 'poisson') {
      distribution(y) <- greta::poisson(lambda = exp(mu))
    } else {
      stop("family must be 'gaussian' or 'poisson'")
    }
  }
  
  # define model
  if (!is.null(z)) {
    greta_model <- greta::model(mu, alpha, beta,
                                sigma_gamma, sigma_main, sigma_bins,
                                ...)
  } else {
    greta_model <- greta::model(mu, alpha, beta,
                                sigma_main, sigma_bins,
                                ...)
  }

  # return model
  list(greta_model = greta_model,
       spline_basis = spline_basis,
       bins = bins)
  
}

# internal function: create greta_fda object
as.greta_fda <- function (model) {
  as_class(model, name = 'greta_fda', type = 'list')
}
