#' A fda_predictor object containing a greta_array for use with function-valued predictor variables
#'
#' @description A \code{fda_predictor} object contains a function regression object
#'   that can be used to model a function-valued predictor with \link[greta]{model}
#' 
#' @rdname fda_predictor
#' 
#' @param x a \code{numeric} vector, \code{matrix}, or \code{data.frame} with a single predictor variable (or a fda_predictor object for is.fda_predictor and print methods)
#' @param bins a \code{integer} or \code{numeric} vector or \code{matrix} of values at which the predictor variable is recorded
#' @param spline_settings a named list of settings to pass to the spline function (see \link[splines]{bs} for details)
#' @param priors a named list of settings for prior distributions (mean, sd, sigma_max)
#' @param object a \code{fda_predictor} object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{fda_predictor}, which can be used in a \link[greta]{model}
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
#' fda_x <- fda_predictor(example_fda_data$x_fun,
#'                              priors = list(
#'                                nu_local = 2,
#'                                nu_global = 2,
#'                                slab_df = 2,
#'                                scale_global = 2,
#'                                slab_scale = 2
#'                              )))
#'                         
#' \dontrun{                 
#' 
#' # define some priors
#' alpha <- normal(mean = 0.0, sd = 1.0, dim = 1)
#' beta <- normal(mean = 0.0, sd = 1.0, dim = 1)
#' sigma_response <- uniform(min = 0.0, max = 5.0, dim = 1)
#' 
#' # define a linear predictor with the fda_predictor object
#' mu <- alpha + beta * example_fda_data$x1 + fda_x$array
#' 
#' # fit a model
#' distribution(example_fda_data$y_scalar) <- normal(mu, sd = sigma_response)
#' 
#' greta_model <- model(mu, alpha, beta, fda_x$beta, sigma_response)
#' }

fda_predictor <- function (x,
                           bins = NULL,
                           spline_settings = list(),
                           priors = list(), ...) {
  
  # test dimensions and class of input
  if (!is.matrix(x)) {
    if (is.data.frame(x)) {
      classes <- apply(x, 2, class)
      if (!all(classes == "numeric")) {
        stop(paste0('the following columns in x are not numeric: ',
                    colnames(x)[which(classes != 'numeric')]))
      }
      x <- as.matrix(x)
      model_type <- 'matrix'
    } else {
      if (is.numeric(x)) {
        model_type <- 'flat'
        if (is.null(bins)) {
          stop('bins must be provided if the predictor variable is flattened')
        }
      } else {
        stop('x must be a numeric matrix, vector, or data.frame')
      }
    }
  } else {
    model_type <- 'matrix'
  }

  # unpack spline settings
  spline_set <- list(basis = "bs",
                     df = 10,
                     degree = 3,
                     intercept = FALSE)
  spline_set[names(spline_settings)] <- spline_settings
  
  # prepare greta_array
  if (model_type == 'matrix') {
    fda_predictor <- build_fda_predictor_matrix(x,
                                                bins,
                                                priors,
                                                spline_set,
                                                ...)
  }
  if (model_type == 'flat') {
    fda_predictor <- build_fda_predictor_flat(x,
                                              bins,
                                              priors,
                                              spline_set,
                                              ...)
  }  
  
  as.fda_predictor(fda_predictor)
  
}

#' @rdname fda_predictor
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is a fda_predictor object
#'   
#' \dontrun{
#' is.fda_predictor(x)
#' }

is.fda_predictor <- function (x) {
  inherits(x, 'fda_predictor')
}

#' @rdname fda_predictor
#'
#' @export
#'
#' @examples
#' 
#' # Print information about a 'fda_predictor' object
#'
#' \dontrun{
#' print(x)
#' }

print.fda_predictor <- function (x, ...) {
  cat(paste0('This is a fda_predictor object\n'))
}

#' @rdname fda_predictor
#'
#' @export
#'
#' @examples
#' 
#' # Summarise a 'fda_predictor' object
#'
#' \dontrun{
#' summary(object)
#' }

summary.fda_predictor <- function (object, ...) {
  cat(paste0('This is a fda_predictor object\n'))
}


# internal function: create greta_array from matrix input data
build_fda_predictor_matrix <- function (x,
                                        bins,
                                        priors,
                                        spline_settings, ...) {
  
  # pull out index counters
  n <- nrow(x)
  nj <- ncol(x)

  # convert x to a greta array
  x <- greta::as_data(x)

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
                                             intercept = spline_settings$intercept,
                                             Boundary.knots = boundary_knots)
  spline_basis <- greta::as_data(spline_basis)

  # setup priors
  prior_set <- list(nu_local = 2,
                    nu_global = 2,
                    slab_df = 2,
                    scale_global = 2,
                    slab_scale = 2)
  prior_set[names(priors)] <- priors

  # setup parameters
  aux1_local <- greta::normal(0, 1, dim = np)
  aux2_local <- greta::inverse_gamma(
    0.5 * prior_set$nu_local,
    0.5 * prior_set$nu_local,
    dim = np
  )
  aux1_global <- greta::normal(0, 1, dim = 1)
  aux2_global <- greta::inverse_gamma(
    0.5 * prior_set$nu_global,
    0.5 * prior_set$nu_global,
    dim = 1
  )
  caux <- greta::inverse_gamma(
    0.5 * prior_set$slab_df,
    0.5 * prior_set$slab_df,
    dim = 1
  )
  lambda <- aux1_local * sqrt(aux2_local)
  tau <- aux1_global * sqrt(aux2_global) * prior_set$scale_global
  cterm <- prior_set$slab_scale * sqrt(caux)
  lambda_tilde <- sqrt((cterm ^ 2 * lambda ^ 2) / 
                         (cterm ^ 2 + tau ^ 2 * lambda ^ 2))
  z <- greta::normal(0, 1, dim = np)
  beta <- z * lambda_tilde * tau

  # define linear predictor
  greta_array <- x %*% (spline_basis %*% beta)

  # return model
  list(array = greta_array,
       beta = beta,
       spline_basis = spline_basis,
       bins = bins)
  
}

# internal function: create greta array from flattened input data
build_fda_predictor_flat <- function (x,
                                      bins,
                                      priors,
                                      spline_settings, ...) {
  
  # pull out index counters
  n <- length(x)

  # convert x and y to greta arrays
  x <- greta::as_data(x)

  # set up spline settings (nspline, nknots, degree)
  boundary_knots <- c(0, max(bins) + 1)
  np <- spline_settings$df
  
  # create spline basis
  spline_basis <- get(spline_settings$basis)(bins,
                                             df = np,
                                             degree = spline_settings$degree,
                                             intercept = spline_settings$intercept,
                                             Boundary.knots = boundary_knots)
  spline_basis <- greta::as_data(spline_basis)

  # setup priors
  prior_set <- list(nu_local = 2,
                    nu_global = 2,
                    slab_df = 2,
                    scale_global = 2,
                    slab_scale = 2)
  prior_set[names(priors)] <- priors

  # setup parameters
  aux1_local <- greta::normal(0, 1, dim = np)
  aux2_local <- greta::inverse_gamma(
    0.5 * prior_set$nu_local,
    0.5 * prior_set$nu_local,
    dim = np
  )
  aux1_global <- greta::normal(0, 1, dim = 1)
  aux2_global <- greta::inverse_gamma(
    0.5 * prior_set$nu_global,
    0.5 * prior_set$nu_global,
    dim = 1
  )
  caux <- greta::inverse_gamma(
    0.5 * prior_set$slab_df,
    0.5 * prior_set$slab_df,
    dim = 1
  )
  lambda <- aux1_local * sqrt(aux2_local)
  tau <- aux1_global * sqrt(aux2_global) * prior_set$scale_global
  cterm <- prior_set$slab_scale * sqrt(caux)
  lambda_tilde <- sqrt((cterm ^ 2 * lambda ^ 2) / 
                         (cterm ^ 2 + tau ^ 2 * lambda ^ 2))
  z <- greta::normal(0, 1, dim = np)
  beta <- z * lambda_tilde * tau
  
  # define linear predictor
  greta_array <- x * (spline_basis %*% beta)

  # return model
  list(array = greta_array,
       beta = beta,
       spline_basis = spline_basis,
       bins = bins)
  
}

# internal function: create a fda_predictor object
as.fda_predictor <- function (x) {
  as_class(x, name = 'fda_predictor', type = 'list')
}
