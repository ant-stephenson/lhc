library(testthat)
library(fs)

source_path <- path_join(c(path_dir(path_dir(getwd())), "R"))
filepath <- path_join(c(source_path, "utility_funcs.r"))

source(filepath)

test_that("pairwise_distance", {
  X0 <- matrix(1:10, nrow=2, ncol=2)
  X1 <- matrix(2:11, nrow=2, ncol=2)
  s <- pairwise_distance(X0, X1)
  expect_equal(s, c(2,2))
})

test_that("inv_model_idx", {
  K <- 5
  midx <- matrix(, nrow=6*K,ncol=2)
    for (j in 1:6) {
      for (k in 1:K) {
        idx <- get_model_idx(j, k, K)
        midx[idx,] <- inv_model_idx(idx, K)
      }
    }
  expected <- matrix(c(rep(1:6, each=5), rep(1:5, 6)), nrow=6*K, ncol=2)  
  expect_equal(midx, expected)
})
