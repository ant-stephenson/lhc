#' Soft margin SVM
#'
#' Fit a soft margin support vector machine for binary classification.
#'
#' @param X An nxd matrix with samples as rows and features as columns.
#' @param y A length-n vector of 0s and 1s indicating true sample classes.
#' @param C Regularisation parameter.
#'
#' @importFrom quadprog solve.QP
#' @export
svm <- function(X, y, C=1){
  #we want to find the values of lambda which minimise (lambda+y)^T * X * X^T * (lambda+y) - <lambda, 1>
  #subject to 0 <= lambda_i <= C 0 and sum lambda_i y_i = 0

  #to use quadprog we need to rewrite the minimisation in the form 0.5 * lambda^T * D * lambda - d^T lambda
  #so to incorporate the y values into X (which we'll call Q)
  n <- nrow(X)
  Q <- sapply(1:n, function(i) y[i]*X[i,])
  D <- crossprod(Q, Q)
  d <- as.matrix(rep(1, n), ncol=1)

  #using a trick to perturb D by a small amount on the diagonal to ensure its positive definite
  D <- D + 1e-6 * diag(n)

  #then we need to represent the constraints in the form A^T * lambda >= b
  #to get each value of lambda > 0 we can set A = I and b = 0
  #to get each value of lambda < C we can add another n rows where A = -I and b =-C
  #to get the sum of lambda_i y_i = 0 we can add a new row to A which is the value of y (the first row)
  #so A is a (2n+1) by n matrix and b is a (2n+1) vector
  #using the argument meq=1 the first constraint is an equality constraint and the rest are inequality
  A <- t(rbind(y, diag(n), -1*diag(n)))
  b <- c(0, rep(0, n), rep(-C, n))

  # call the QP solver:
  sol <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec=b, meq=1)
  lambda <- sol$solution

  # convert lambda back to w and w0
  w <- rep(0, ncol(X))
  for(i in 1:n){
    w <- w + lambda[i] * y[i] %*% X[i,]
  }

  # support vectors have lambda != 0
  svs <- which(lambda > 1e-4)

  # find w_0 using the support vectors
  w0 <- 0
  for(i in svs){
    w0 <- w0 + y[i] - w %*% X[i,]
  }
  w0 <- w0 / length(svs)

  # return a function to predict class of new points
  f <- function(X_test){
    y <- X_test %*% t(w) + rep(w0, nrow(X_test))
    return(sign(y))
  }
}

## Test simple SVM on 2D dataset
# #Generate simple dataset
# n <- 500
# y <- c(rep(1, n/2), rep(-1, n/2))
# x_1 <- c(rnorm(n/2, 2, 2), rnorm(n/2, 8, 3))
# x_2 <- c(rnorm(n/2, 3, 2), rnorm(n/2, 10, 3))
# X <- matrix(c(x_1, x_2), ncol=2)
#
# #shuffle
# ind <- sample(1:n, n)
# y <- y[ind]
# X <- X[ind, ]
# plot(X[,1], X[,2], col=as.factor(y), pch=16)
#
# #plot coloured by predicted class
# model <- svm(X, y, C=1)
# y_pred <- model(X)
# plot(X[,1], X[,2], col=as.factor(y_pred), pch=16)



#' Kernel SVM
#'
#' Fit a kernel support vector machine for binary classification.
#'
#' @param X An nxd matrix with samples as rows and features as columns.
#' @param y A length-n vector of 0s and 1s indicating true sample classes.
#' @param C Regularisation parameter.
#' @param ckernel Kernel function to be used.
#' @param ... Additional parameters passed to ckernel.
#'
#' @importFrom quadprog solve.QP
#' @export
kernel_svm <- function(X, y, C, ckernel, ...){
  #we want to find the values of lambda which minimise (lambda+y)^T * K * (lambda+y) - <lambda, 1>
  #subject to lambda_i >= 0 and sum lambda_i y_i = 0

  #to use quadprog we need to rewrite the minimisation in the form 0.5 * lambda^T * D * lambda - d^T lambda
  #so D= yy^T K
  n <- nrow(X)
  K <- calc_K(X, ckernel, ...)
  D <- outer(y, y) * K
  #using a trick to perturb D by a small amount on the diagonal to ensure its positive definite
  D <- D + 1e-5 * diag(n)
  #get the vector 1s for d
  d <- as.matrix(rep(1, n), ncol=1)

  #then we need to represent the constraints in the form A^T * lambda >= b
  #to get each value of lambda > 0 we can set A = I and b = 0
  #to get each value of lambda < C we can add another n rows where A = -I and b =-C
  #to get the sum of lambda_i y_i = 0 we can add a new row to A which is the value of y (the first row)
  #so A is a (2n+1) by n matrix and b is a (2n+1) vector
  #using the argument meq=1 the first constraint is an equality constraint and the rest are inequality
  A <- t(rbind(y, diag(n), -1*diag(n)))
  b <- c(0, rep(0, n), rep(-C, n))

  # call the QP solver:
  sol <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec=b, meq=1)
  lambda <- sol$solution

  #find support vectors
  svs <- which(lambda > 1e-4)

  #could have a check if too many support vectors (likely means overfitting and predictions will be slow)

  #the equivalent of w^T x_i for kernel svm is sum_j lambda[j]*y[j]*k(x_i, x_j)
  #but for non-support vectors lambda = 0, so really we only need to sum over the support vec
  #define a function to compute this for value for a sample
  val <- function(x_i){
    v <- 0
    for(j in svs){
      v <- v + lambda[j] * y[j] * ckernel(x_i, X[j,], ...)
    }
    return(v)
  }

  #choose any support vector and find w0 using equation 7.17 from PRML
  #(equation 7.18 says its better to take an average across all svs, but this is more computationally intensive)
  # s <- which.max(lambda)
  # w0 <- y[s] - val(X[s,])
  w0 <- 0
  for(s in svs){
    w0 <- w0 + y[s] - val(X[s,])
  }
  w0 <- w0 / length(svs)

  # return a function to predict class of new points
  f <- function(X_test){
    y <- apply(X_test, 1, function(x_i) val(x_i) + w0)
    return(sign(y))
  }
}


## Test the kernel svm on a non linear dataset
# #Generate data
# n <- 500
# c1_x1 <- rnorm(n/2, 0, 1.2)
# c1_x2 <- rnorm(n/2, 0, 1.2)
#
# c2_x1 <- seq(-4, 4, length.out = n/4)
# c2_x2 <- sqrt(16 - c2_x1^2)
# c2_x1 <- c(c2_x1, c2_x1) + rnorm(n/2, 0, 0.8)
# c2_x2 <- c(c2_x2, -c2_x2) + rnorm(n/2, 0, 0.8)
#
# X <- matrix(c(c1_x1, c2_x1, c1_x2, c2_x2), nrow=n)
# y <- c(rep(-1, n/2), rep(1, n/2))
#
# #shuffle
# ind <- sample(1:n, n)
# y <- y[ind]
# X <- X[ind, ]
# plot(X[,1], X[,2], col=as.factor(y), pch=16) #true class
#
# #polynomial kernel
# model <- kernel_svm(X, y, C=1, ckernel=poly_kernel, b=3)
# y_pred <- model(X)
# plot(X[,1], X[,2], col=as.factor(y_pred), pch=16) #predicted class
# #RBF kernel
# sigma <- avg_median_pairwise_distance(X)
# model <- kernel_svm(X, y, C=1, ckernel=rbf_kernel, sigma=sigma)
# y_pred <- model(X)
# plot(X[,1], X[,2], col=as.factor(y_pred), pch=16)
# #linear kernel
# model <- kernel_svm(X, y, C=1, ckernel=lin_kernel)
# y_pred <- model(X)
# plot(X[,1], X[,2], col=as.factor(y_pred), pch=16)
# #trigonometric kernel
# model <- kernel_svm(X, y, C=1, ckernel=trig_kernel, b=6)
# y_pred <- model(X)
# plot(X[,1], X[,2], col=as.factor(y_pred), pch=16)
