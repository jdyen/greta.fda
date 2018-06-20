#' Example data sets for functional data analysis
#'
#' Example data set with a functional response variable, two covariates, and one grouping variable
#'
#' @rdname example_datasets
#'
#' @docType data
#'
#' @usage example_fda_data
#'
#' @format A list with three elements: response matrix (y), covariates (x1, x2), and a grouping variable (z1)
#'
#' @keywords datasets
#'
#' @examples
#' library(greta.fda)
#' names(example_fda_data)
#' head(example_fda_data$y)
'example_fda_data'