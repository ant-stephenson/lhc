#' Soft margin SVM
#'
#' Fit a soft margin support vector machine for binary classification.
#'
#' @param X An nxd matrix with samples as rows and features as columns.
#' @param y A length-n vector of -1s and 1s indicating true sample classes.
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
  Q <- tcrossprod(X, X)
  D <- outer(y, y) * Q
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

  # support vectors have lambda != 0
  svs <- which(lambda > 0.01)

  # convert lambda back to w and w0
  w <- rep(0, ncol(X))
  for(i in 1:n){
    w <- w + lambda[i] * y[i] %*% X[i,]
  }

  w0 <- 0
  for(i in svs){
    w0 <- w0 + y[i] - w %*% X[i,]
  }
  w0 <- w0 / length(svs)

  # return a function to predict class of new points
  f <- function(X_test){
    y <- X_test %*% t(w) + rep(w0, nrow(X_test))
    return(as.numeric(sign(y)))
  }
}


#' Kernel SVM
#'
#' Fit a kernel support vector machine for binary classification.
#'
#' @param X An nxd matrix with samples as rows and features as columns.
#' @param y A length-n vector of -1s and 1s indicating true sample classes.
#' @param C Regularisation parameter.
#' @param ckernel Kernel function with hyperparamters set
#'
#' @importFrom quadprog solve.QP
#' @export
kernel_svm <- function(X, y, C, ckernel){
  #we want to find the values of lambda which minimise (lambda+y)^T * K * (lambda+y) - <lambda, 1>
  #subject to lambda_i >= 0 and sum lambda_i y_i = 0

  #to use quadprog we need to rewrite the minimisation in the form 0.5 * lambda^T * D * lambda - d^T lambda
  #so D= yy^T K
  n <- nrow(X)
  K <- calc_K(X, ckernel)
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
  svs <- which(lambda > 0.01)

  #could have a check if too many support vectors (likely means overfitting and predictions will be slow)

  #the equivalent of w^T x_i for kernel svm is sum_j lambda[j]*y[j]*k(x_i, x_j)
  #but for non-support vectors lambda = 0, so really we only need to sum over the support vec
  #define a function to compute this for value for a sample
  val <- function(x_i){
    v <- 0
    for(j in svs){
      v <- v + lambda[j] * y[j] * ckernel(x_i, X[j,])
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

