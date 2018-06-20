context('fda_response class')

test_that('fda_response formula works', {
  
  expect_ok(fda_response(y ~ x1 + x2 + (1 | z1), data = example_fda_data))
  expect_ok(fda_response(y ~ x1 + x2, data = example_fda_data)) 
  
  data_set <- example_fda_data
  data_set$x1 <- rep(data_set$x1, times = ncol(data_set$y))
  data_set$x2 <- rep(data_set$x2, times = ncol(data_set$y))
  data_set$z1 <- rep(data_set$z1, times = ncol(data_set$y))
  data_set$y <- c(data_set$y)
  expect_ok(fda_response(y ~ x1 + x2 + (1 | z1), data = data_set,
                         bins = rep(seq_len(ncol(example_fda_data$y)), each = nrow(example_fda_data$y))))
  expect_ok(fda_response(y ~ x1 + x2, data = data_set,
                         bins = rep(seq_len(ncol(example_fda_data$y)), each = nrow(example_fda_data$y))))
   
}) 
