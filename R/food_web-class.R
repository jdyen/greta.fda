#' Create a food_web object to use in a trophic projection
#'
#' @description A 'food web' object stores the trophic dynamics of a set of potentially interacting species.
#' It is a sub-component of a \link[trophic]{trophic_dynamics} object
#' 
#' @rdname food_web
#' 
#' @param interaction_matrix a matrix identifying (weighted) trophic links between species pairs (rows eat columns)
#' @param x an object to print or test as a food_web object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{food_web}
#' 
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix components induced.subgraph shortest.paths V E
#' @importFrom graphics plot
#' 
#' @examples
#' 
#' library(trophic)
#' 
#' # Construct the food_web object
#' 
#' test_fw <- build_food_web(interaction_matrix = food_web)

build_food_web <- function (interaction_matrix, ...) {
  
  # check food web for cannibalism
  if (!all(almost_equal(diag(interaction_matrix), 0))) {
    stop("A diagonal element of the interaction_matrix is not zero but cannibalism is not supported")
  }
  
  # check if a symmetric matrix has been provided
  if (isSymmetric(interaction_matrix)) {
    stop("The interaction_matrix is symmetric but must not contain loops")
  }
  
  # check food web for loops
  upper_matrix <- interaction_matrix * upper.tri(interaction_matrix)
  lower_matrix <- interaction_matrix * lower.tri(interaction_matrix)
  if (!all(almost_equal(upper_matrix, 0)) & !all(almost_equal(lower_matrix, 0))) {
    stop("There are non-zero entries in the lower and upper diagonal of the interaction_matrix but this matrix should be sorted into lower or upper diagonal form")
  }
  
  # convert food web to lower triangular form
  if (!all(almost_equal(upper_matrix, 0))) {
    interaction_matrix <- upper_matrix[rev(seq_len(nrow(interaction_matrix))), rev(seq_len(ncol(interaction_matrix)))]
  } else {
    interaction_matrix <- lower_matrix
  }

  # add column and row names
  if (is.null(colnames(interaction_matrix)) & is.null(rownames(interaction_matrix))) {
    colnames(interaction_matrix) <- letters[seq_len(ncol(interaction_matrix))]
    rownames(interaction_matrix) <- colnames(interaction_matrix)
  } else {
    if (is.null(rownames(interaction_matrix))) {
      rownames(interaction_matrix) <- colnames(interaction_matrix)
    } else {
      colnames(interaction_matrix) <- rownames(interaction_matrix)
    }
  }
  
  # identify if food web is fixed or stochastic
  type <- "stochastic"
  if (all(c(interaction_matrix) %in% c(0, 1))) {
    type <- "fixed"
  }
  
  # create food_web object
  food_web <- list(interaction_matrix = interaction_matrix,
                   nsp = ncol(interaction_matrix),
                   type = type)

  # return food_web object with class definition
  as.food_web(food_web)
  
}

#' @rdname food_web
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'food_web'
#'   
#' \dontrun{
#' is.food_web(x)
#' }

is.food_web <- function (x) {
  inherits(x, 'food_web')
}

#' @rdname food_web
#'
#' @export
#'
#' @examples
#' 
#' # Print information about the 'food_web' object
#'
#' \dontrun{
#' print(x)
#' }

print.food_web <- function (x, ...) {
  cat(paste0("This is a ", x$type, " food_web object with ", x$nsp, " species"))
}

#' @rdname food_web
#'
#' @export
#'
#' @examples
#' 
#' # Plot a 'food_web' object
#'
#' \dontrun{
#' plot(x)
#' }

plot.food_web <- function (x, ...) {

  # create an igraph object from the adjacency matrix
  fw_graph <- igraph::graph_from_adjacency_matrix(t(x$interaction_matrix),
                                                  weighted = TRUE,
                                                  mode = "directed")
  
  # labels need to ignore disconnected nodes
  label_subgraph <- rownames(x$interaction_matrix)
  label_subgraph <- label_subgraph[!((apply(x$interaction_matrix, 1, sum) == 0) &
                                       (apply(x$interaction_matrix, 2, sum) == 0))]
  
  # plot the adjacency matrix
  plot(fw_graph,
       layout = layout_fw,
       vertex.label = label_subgraph,
       ...)
  
}

# internal function: create food_web object
as.food_web <- function (food_web) {
  as_class(food_web, name = "food_web", type = "list")
}
