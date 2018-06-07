# test functions

# simulate a simple food web without stochasticity to test calculation
fixed_projection <- function(trophic_dynamics, primary_producers) {
  
  # only run for a single food web model
  i <- 1
  nsim <- 2
  
  # unpack indices and food webs
  nsp <- trophic_dynamics$food_web[[i]]$nsp
  if (trophic_dynamics$ntrophic > 1) {
    food_web <- trophic_dynamics$food_web[[i]]
  } else {
    food_web <- trophic_dynamics$food_web[[1]]
  }
  if (trophic_dynamics$nefficiency > 1) {
    efficiency_matrix <- trophic_dynamics$efficiency_matrix[[i]]
  } else {
    efficiency_matrix <- trophic_dynamics$efficiency_matrix[[1]]
  }
  if (trophic_dynamics$ndominance > 1) {
    dominance_matrix <- trophic_dynamics$dominance_matrix[[i]]
  } else {
    dominance_matrix <- trophic_dynamics$dominance_matrix[[1]]
  }
  
  # match primary producer names to food web
  spnames <- rownames(food_web$interaction_matrix)
  primary_producer_id <- match(primary_producers$names, spnames)
  if (any(is.na(primary_producer_id))) {
    primary_producer_id <- seq_len(length(primary_producers))
    warning(paste0("names of primary_producers do not match nodes in food_web; primary_producers are assumed to be the first ",
                   length(primary_producers),
                   " nodes of food_web"))
  }
  
  # create production output matrix
  production <- matrix(0, nrow = nsp, ncol = nsim)
  rownames(production) <- spnames
  
  # initialise production matrix
  production[primary_producer_id, ] <- matrix(rep(unlist(primary_producers$mean), times = nsim),
                                              nrow = primary_producers$n)
  
  # calculate dominance weights
  n_split <- apply(dominance_matrix$dominance, 2, sum)
  if (any(n_split == 0)) {
    n_split[n_split == 0] <- 1
  }
  dominance_weights <- sweep(dominance_matrix$dominance,
                             2,
                             n_split,
                             "/")
  dominance_weights <- as.matrix(dominance_weights)
  
  # calculate efficiency over all iterations
  efficiency <- array(rep(efficiency_matrix$mean, times = nsim),
                      dim = c(nsp, nsp, nsim))
  
  # calculate amount of biomass in each node (working from primary consumers up the food web)
  nodeorder <- seq_len(food_web$nsp)
  nodeorder <- nodeorder[apply(food_web$interaction_matrix, 1, sum) > 0]
  for (node in nodeorder) {
    for (sim in seq_len(nsim)) {
      for (sp in seq_len(nsp)) {
        production[node, sim] <- production[node, sim] +
          dominance_weights[node, sp] * (efficiency[node, sp, sim] * production[sp, sim])
      }
    }
  }
  
  
  production
  
}
