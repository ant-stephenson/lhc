#' Define polynomial kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @param b order of polynomial
#' @return scalar of (1+x^Tx)^b
#' @export
poly_kernel <- function(x_i, x_j, b) {
  return((t(x_i) %*% x_j + 1)^b)
}

#' Define trigonometric kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @param b order of polynomial
#' @export
trig_kernel <- function(x_i, x_j, b=0) {
  output <- 0
  for (k in 0:b) {
    output <- output + cos(k * (x_i - x_j))
  }
  return(sum(output))
}

#' Define linear kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @export
lin_kernel <- function(x_i, x_j) {
  return(t(x_i) %*% x_j)
}

#' RBF kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @param sigma bandwidth hyperparameter
#' @export
rbf_kernel <- function(x_i, x_j, sigma) {
  return(exp(-t(x_i - x_j) %*% (x_i - x_j) / (2 * sigma^2)))
}

#' Function factory to partially call kernel function to return the kernel function with it's hyperparameters set
#' @param ckernel kernel function with args (x_i,x_j,hyper)
#' @return function with args (x_i,x_j)
#' @export
tuned_kernel <- function(ckernel, ...) {
  function(x_i, x_j){
    ckernel(x_i, x_j, ...)
  }
}

#' Compute the kernel matrix over the (training) set X
#' @param X covariate matrix (nxd)
#' @param ckernel kernel function with args (x_i, x_j)
#' @return K kernel matrix (nxn)
#' @export
calc_K <- function(X, ckernel){
  n <- nrow(X)
  K <- matrix(NA, nrow=n, ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      K[i,j] <- ckernel(X[i,], X[j,])
    }
  }
  return(K)
}
