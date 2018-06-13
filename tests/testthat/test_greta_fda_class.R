context('greta_fda class')

test_that('greta_fda formula works', {
  
  expect_ok(greta_fda(y ~ x1 + x2 + (1 | z1), data = greta_fda_data,
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1)))
  expect_ok(greta_fda(y ~ x1 + x2, data = greta_fda_data,
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1))) 
  
  data_set <- greta_fda_data
  data_set$y <- round(data_set$y) + 50
  expect_ok(greta_fda(y ~ x1 + x2 + (1 | z1), data = data_set,
                      family = 'poisson',
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1))) 
  expect_ok(greta_fda(y ~ x1 + x2, data = data_set,
                      family = 'poisson', 
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1)))
  
  data_set <- greta_fda_data
  data_set$x1 <- rep(data_set$x1, times = ncol(data_set$y))
  data_set$x2 <- rep(data_set$x2, times = ncol(data_set$y))
  data_set$z1 <- rep(data_set$z1, times = ncol(data_set$y))
  data_set$y <- c(data_set$y)
  expect_ok(greta_fda(y ~ x1 + x2 + (1 | z1), data = data_set,
                      family = 'gaussian',
                      bins = rep(seq_len(ncol(greta_fda_data$y)), each = nrow(greta_fda_data$y)),
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1)))
  expect_ok(greta_fda(y ~ x1 + x2, data = data_set,
                      family = 'gaussian', 
                      bins = rep(seq_len(ncol(greta_fda_data$y)), each = nrow(greta_fda_data$y)),
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1)))

  data_set$y <- round(data_set$y) + 50
  expect_ok(greta_fda(y ~ x1 + x2 + (1 | z1), data = data_set,
                      family = 'poisson', 
                      bins = rep(seq_len(ncol(greta_fda_data$y)), each = nrow(greta_fda_data$y)),
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1)))
  expect_ok(greta_fda(y ~ x1 + x2, data = data_set,
                      family = 'poisson', 
                      bins = rep(seq_len(ncol(greta_fda_data$y)), each = nrow(greta_fda_data$y)),
                      greta_settings = list(n_samples = 100, warmup = 100, chains = 1)))
  
}) 
