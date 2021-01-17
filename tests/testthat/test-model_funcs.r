library(testthat)
library(fs)

source_path <- path_join(c(path_dir(path_dir(getwd())), "R"))
filepath <- path_join(c(source_path, "lr_funcs.r"))

source(filepath)

test_logistic_reg_implementation <- function(X ,y) {
    glm_b <- glm(y~. - 1, family="binomial", data = X)
    b_hat <- logistic_reg(as.matrix(X), y, lambda=0)

    # check how similar
    max_diff <- max(abs(coef(glm_b) - b_hat))
    return(max_diff)
    stopifnot(max_diff < 1e-6)
}

test_that("logistic regression", {
  X <- model.matrix(~1 + x1 + x2, data.frame(x1=rnorm(1000), x2=rnorm(1000)))
  b <- rnorm(3)
  y <- rbinom(1000, 1, logisticf(X %*% b))
  max_diff <- test_logistic_reg_implementation(data.frame(X), y)
  expect_lt(max_diff, 1e-6)
})
