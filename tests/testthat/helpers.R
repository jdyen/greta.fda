# test functions
expect_ok <- function (expr)
  expect_error(expr, NA)