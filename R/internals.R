# internal helper functions for greta.fda R package

# set an object class
as_class <- function (object, name, type = c("function", "list")) {
  
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  
  object
  
}

# test if values are equal to zero with a non-zero tolerance
almost_equal <- function(x, y, tolerance = 1e-10) {
  
  diff <- abs(x - y)
  mag <- pmax(abs(x), abs(y))
  
  ifelse(mag > tolerance, (diff / mag) <= tolerance, diff <= tolerance)
  
}

# test if multiple values are all equal
all_equal <- function (..., tolerance = 1e-10) {
  
  test <- list(...)
  
  if (length(test) > 1) {
    out <- rep(NA, length(test))
    for (i in seq_along(test)) {
      out[i] <- almost_equal(test[[i]], test[[1]], tolerance = tolerance)
    }
  } else {
    out <- TRUE
  }
  
  all(out)
  
}

# define shrinkage prior
define_shrinkage_prior <- function(priors, np) {
  
  # setup parameters
  aux1_local <- greta::normal(0, 1, dim = np)
  aux2_local <- greta::inverse_gamma(
    0.5 * priors$nu_local,
    0.5 * priors$nu_local,
    dim = np
  )
  aux1_global <- greta::normal(0, 1, dim = 1)
  aux2_global <- greta::inverse_gamma(
    0.5 * priors$nu_global,
    0.5 * priors$nu_global,
    dim = 1
  )
  caux <- greta::inverse_gamma(
    0.5 * priors$slab_df,
    0.5 * priors$slab_df,
    dim = 1
  )
  lambda <- aux1_local * sqrt(aux2_local)
  tau <- aux1_global * sqrt(aux2_global) * priors$scale_global
  cterm <- priors$slab_scale * sqrt(caux)
  lambda_tilde <- sqrt((cterm ^ 2 * lambda ^ 2) / 
                         (cterm ^ 2 + tau ^ 2 * lambda ^ 2))
  z <- greta::normal(0, 1, dim = np)
  z * lambda_tilde * tau
  
}
