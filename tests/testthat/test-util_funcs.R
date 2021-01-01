library(testthat)
source("../../utility_funcs.r")

test_that("pairwise distance", {
  X0 <- matrix(1:10, nrow=2, ncol=2)
  X1 <- matrix(2:11, nrow=2, ncol=2)
  s <- pairwise_distance(X0, X1)
  expect_equal(s, c(5,5))
})
