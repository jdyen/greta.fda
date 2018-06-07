context('food_web class')

test_that('build_food_web works', {
  
  test_fw <- matrix(rbinom(16, size = 1, p = 0.5), ncol = 4)
  test_fw <- test_fw * lower.tri(test_fw)
  
  # upper triangular
  expect_silent(build_food_web(t(test_fw)))
  
  # lower triangular
  expect_silent(build_food_web(test_fw))
  
  # with colnames
  colnames(test_fw) <- letters[1:4]
  expect_silent(build_food_web(test_fw))
  
  # with rownames
  colnames(test_fw) <- NULL
  rownames(test_fw) <- letters[1:4]
  expect_silent(build_food_web(test_fw))
  
  # with colnames and rownames
  colnames(test_fw) <- letters[1:4]
  expect_silent(build_food_web(test_fw))

  # continuous food web
  test_fw <- test_fw * matrix(runif(16, min = 0.5, max = 1.5), ncol = 4)
  expect_silent(build_food_web(test_fw))
  
})

test_that('print works', {

  test_fw <- matrix(rbinom(16, size = 1, p = 0.5), ncol = 4)
  test_fw <- test_fw * lower.tri(test_fw)
  test_fw <- build_food_web(test_fw)
  
  # print method
  expected_output <- "This is a fixed food_web object with 4 species"
  expect_output(print(test_fw), expected_output)

})

test_that('is.food_web works', {

  test_fw <- matrix(rbinom(16, size = 1, p = 0.5), ncol = 4)
  test_fw <- test_fw * lower.tri(test_fw)
  test_fw <- build_food_web(test_fw)
  
  # is.food_web expect TRUE
  expect_true(is.food_web(test_fw))

  # is.food_web expect FALSE
  expect_false(is.food_web(rnorm(100)))
  
  # errors with no argument
  expect_error(is.food_web())
  
})

test_that('plot works', {
  
  # plot works for any food_web object
  test_fw <- matrix(rbinom(16, size = 1, p = 0.5), ncol = 4)
  test_fw <- test_fw * lower.tri(test_fw)
  test_fw <- build_food_web(test_fw)
  expect_silent(plot(test_fw))

})

test_that('food_web errors are correct', {
  
  source('helpers.R')

  test_fw <- matrix(rbinom(16, size = 1, p = 0.5), ncol = 4)
  test_fw <- test_fw * lower.tri(test_fw)

  # errors if diag(x) > 0
  test_fw[2, 2] <- 0.5
  expect_error(build_food_web(test_fw))
  
  # errors if x is symmetric
  test_fw[2, 2] <- 0.0
  test_fw[upper.tri(test_fw)] <- t(test_fw)[upper.tri(test_fw)]
  expect_error(build_food_web(test_fw))
  
  # errors if x has loops (not symmetric)
  test_fw <- test_fw * lower.tri(test_fw)
  test_fw[1, 3] <- 0.5
  expect_error(build_food_web(test_fw))
  
})


