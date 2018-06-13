#' A greta_fda regression model with function-valued data
#'
#' @description A \code{greta_fda} object contains a function regression model fitted with \code{greta_fda}
#' 
#' @rdname greta_fda
#' 
#' @param formula formula describing the model to be fitted, in the format used by \link[lme4]{lmer}
#' @param data a named list containing the variables in \code{formula}
#' @param family a GLM family passed as a \code{character}, see \link[stats]{family} (currently only gaussian (default) and poisson families are implemented)
#' @param link link function to be used for model family, see \link[stats]{family} for details (currently not implemented)
#' @param y a \code{matrix} or \code{data.frame} with the response variable (one row per observation)
#' @param x a \code{matrix} or \code{data.frame} of predictor variables
#' @param z a \code{matrix} or \code{data.frame} of random effects variables
#' @param bins a vector of values at which the response variable is recorded
#' @param greta_settings a named list of values to pass to the inference method (see \link[greta]{inference} for details)
#' @param spline_settings a named list of settings to pass to the spline function (see \link[splines]{bs} for details)
#' @param priors a named list of prior distributions in the format of \link[greta]{distributions}
#' @param errors a character denoting the type of errors in a matrix model; currently only 'iid' and 'ar1' are implemented
#' @param model a fitted \code{greta_fda} model
#' @param type for fitted and predict; link or response scale?
#' @param probs for fitted; which quantiles to calculate?
#' @param newdata a list of data for which posterior predictions should be generated
#' @param re.form form of random effects in posterior predictions (not implemented)
#' @param fun function to apply to posterior predictions
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{greta_fda}, which has associated `print`, `plot`, and `summary` methods
#' 
#' @export
#' 
#' @import greta
#' @importFrom splines bs
#' @importFrom graphics plot
#' @importFrom stats terms.formula delete.response fitted quantile as.formula coef model.matrix
#' @importFrom abind abind
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
                             errors = errors,
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
                               family = 'gaussian',
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

  # unpack greta_settings
  greta_set <- list(sampler = greta::hmc(),
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
  if (model_type == 'matrix') {
    greta_model <- build_greta_fda_matrix(y, x, z,
                                          family,
                                          link,
                                          bins,
                                          priors,
                                          errors,
                                          spline_set,
                                          ...)
  }
  if (model_type == 'flat') {
    greta_model <- build_greta_fda_flat(y, x, z,
                                        family,
                                        link,
                                        bins,
                                        priors,
                                        errors,
                                        spline_set,
                                        ...)
  }

  # add initial values if not specified
  greta_set$initial_values <- rep(0.0, length(greta_model$greta_model$dag$example_parameters()))
  
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
  bayes_r2 <- NULL

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

coef.greta_fda <- function (x, 
                            ...,
                            type = c('link', 'response'),
                            probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)) {
  
  # extract coefficients
  param_mean <- lapply(x$samples, function(z) apply(z, 2, mean))
  param_quants <- lapply(x$samples, function(z) t(apply(z, 2, quantile, p = probs)))
  alpha_mean <- lapply(param_mean, function(z) matrix(z[grep('alpha', names(z))],
                                                      nrow = 1))
  alpha_mean <- abind::abind(alpha_mean, along = 3)
  alpha_quants <- lapply(param_quants, function(z) z[grep('alpha', rownames(z)), ])
  alpha_quants <- abind::abind(alpha_quants, along = 3)
  beta_mean <- lapply(param_mean, function(z) matrix(z[grep('beta', names(z))],
                                                     nrow = ncol(x$data$x)))  
  beta_mean <- abind::abind(beta_mean, along = 3)
  beta_quants <- lapply(param_quants, function(z) array(z[grep('beta', rownames(z)), ],
                                                         dim = c(ncol(x$data$x),
                                                                 nrow(x$spline_basis),
                                                                 length(probs))))
  beta_quants <- abind::abind(beta_quants, along = 4)
  alpha <- apply(alpha_mean, c(1, 2), mean)
  beta <- apply(beta_mean, c(1, 2), mean)
  alpha_quants <- apply(alpha_quants, c(1, 2), mean)
  beta_quants <- apply(beta_quants, c(1, 2, 3), mean)
  alpha_link <- alpha %*% as.matrix(x$spline_basis)
  beta_link <- beta %*% as.matrix(x$spline_basis)
  alpha_quants_link <- t(as.matrix(x$spline_basis)) %*% alpha_quants
  beta_quants_link <- array(NA, dim = c(ncol(x$spline_basis),
                                        ncol(x$data$x),
                                        length(probs)))
  for (i in seq_len(dim(beta_quants)[3])) {
    beta_quants_link[, , i] <- t(as.matrix(x$spline_basis)) %*% t(beta_quants[, , i])
  }
  
  # return some summary of this
  list(alpha_mean = alpha_link,
       beta_mean = beta_link,
       alpha_quants = alpha_quants_link,
       beta_quants = beta_quants_link,
       bins = x$data$bins)

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

fitted.greta_fda <- function (x, 
                              ...,
                              type = c('link', 'response'),
                              probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)) {
  
  # calculate fitted values
  param_mean <- lapply(x$samples, function(z) apply(z, 2, mean))
  param_quants <- lapply(x$samples, function(z) t(apply(z, 2, quantile, p = probs)))
  fitted_mean <- lapply(param_mean,
                        function(z) z[grep('mu', names(param_mean[[1]]))])
  fitted_mean <- do.call('rbind', fitted_mean)
  fitted_mean <- apply(fitted_mean, 2, mean)
  fitted_quants <- lapply(param_quants,
                        function(z) z[grep('mu', rownames(param_quants[[1]])), ])
  fitted_quants <- abind::abind(fitted_quants, along = 3)
  fitted_quants <- apply(fitted_quants, c(1, 2), mean)
  
  # work out class of response variable
  if (is.matrix(x$data$y)) {
    dim_set <- ncol(x$data$y)
  } else {
    dim_set <- 1
  }
  
  # reformat mean to same dims as response
  fitted_mean <- matrix(fitted, ncol = dim_set)
  
  # return outputs
  list(mean = fitted_mean,
       quantiles = fitted_quants)
  
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

predict.greta_fda <- function (x, ..., newdata = NULL, type = c('link', 'response'),
                               re.form = NULL, fun = NULL) {
  
  # fill data (check lme4 method for this)
  if (is.null(newdata)) {
    newdata <- x$data
  }
  
  # calculate coefficient
  coefs <- coef(x)
  
  # predict outputs
  out <- NULL
  #out <- coefs$alpha + newdata$x %*% coefs$beta
  
  # add random effects based on re.form
  
  # set link function
  
  
  # return outputs
  out
  
}


# internal function: create greta model from matrix input data
build_greta_fda_matrix <- function (y, x, z,
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
  
  # convert x and y to greta arrays
  x <- greta::as_data(x)
  y <- greta::as_data(y)
  
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
  
  # setup likelihood
  if (!(family %in% c('gaussian', 'poisson', 'binomial'))) {
    stop("family must be 'gaussian' or 'poisson'")
  }
  if (family == 'gaussian') {
    greta::distribution(y) <- greta::normal(mean = mu, sd = sigma_main)
  }
  if (family == 'poisson') {
    greta::distribution(y) <- greta::poisson(lambda = exp(mu))
  }  
  if (family == 'binomial') {
    greta::distribution(y) <- greta::binomial(size = 1,
                                              prob = greta::icloglog(mu))
  } 
  
  # define model
  if (!is.null(z)) {
    gamma_vec <- do.call('c', gamma)
    greta_model <- greta::model(mu,
                                alpha, beta, gamma_vec,
                                sigma_gamma, sigma_main, sigma_bins)
  } else {
    greta_model <- greta::model(mu,
                                alpha, beta,
                                sigma_main, sigma_bins)
  }

  # return model
  list(greta_model = greta_model,
       spline_basis = spline_basis,
       bins = bins)
  
}

# internal function: create greta model from flattened input data
build_greta_fda_flat <- function (y, x, z,
                                  family,
                                  link,
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
  sigma_main <- greta::uniform(min = 0.0, max = 5.0, dim = 1)
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
  mu <- t(alpha %*% spline_basis) + rowSums(x * t(beta %*% spline_basis))
  if (!is.null(z)) {
    for (rand in seq_len(nt)) {
      mu <- mu + rowSums(gamma[[rand]][z[, rand], ] * t(spline_basis))
    }
  }
  
  # setup likelihood
  if (!(family %in% c('gaussian', 'poisson', 'binomial'))) {
    stop("family must be 'gaussian' or 'poisson'")
  }
  if (family == 'gaussian') {
    greta::distribution(y) <- greta::normal(mean = mu, sd = sigma_main)
  }
  if (family == 'poisson') {
    greta::distribution(y) <- greta::poisson(lambda = exp(mu))
  }  
  if (family == 'binomial') {
    greta::distribution(y) <- greta::binomial(size = 1,
                                              prob = greta::icloglog(mu))
  } 
  
  # define model
  if (!is.null(z)) {
    gamma_vec <- do.call('c', gamma)
    greta_model <- greta::model(mu,
                                alpha, beta, gamma_vec,
                                sigma_gamma, sigma_main,
                                ...)
  } else {
    greta_model <- greta::model(mu,
                                alpha, beta,
                                sigma_main,
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
