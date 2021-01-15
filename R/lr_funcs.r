source("utility_funcs.R")

#' Compute newton step
#' @param grad gradient vector
#' @param H Hessian matrix
#' @return step to take
newton_step <- function(grad, H) {
  H_svd <- svd(H)
  H_inv <- H_svd$v %*% diag(1/H_svd$d) %*% t(H_svd$u)

  deltab <- -H_inv %*% grad
  return(deltab)
}
#' backtracking linesearch to find "optimal" step size
#' @param f function we're minimising
#' @param gradf gradient of f
#' @param x parameter we're optimising over
#' @param deltax newton step
#' @param alpha linesearch parameter
#' @param beta linesearch update parameter
backtrack_linesearch <- function(f, gradf, x, deltax, alpha, beta) {
  linesearch_cond <- function(step) {
    f(x + step*deltax) > f(x) + alpha * step * t(gradf(x)) %*% deltax 
  }
  max_iter <- 10
  iter <- 0
  step <- 1
  while (linesearch_cond(step) & iter < max_iter) {
    step <- beta * step
    iter <- iter + 1
  }
  return(step)
}

#' Fit a logistic regression model by IRWLS - same thing glm does
#'
#' @param X covariate matrix
#' @param y response vector
#' @param r [Optional] weight vector
#' @param lambda [Optional] L2 regularisation parameter
#' @return b vector of coefficients
library(Matrix)
logistic_reg <- function(X, y, lambda = 0) {
    invlink <- logistic
    dinvlink <- function(x) exp(-x)/(1+exp(-x))^2
    loss <- function(b) sum(log(1 + exp((-1)^y*(X %*% b))))
    gradloss <- function(b) -t(X) %*% (y - invlink(X %*% b)) + 2 * lambda * b

    d <- ncol(X)
    n <- nrow(X)

    b <- matrix(0, d, 1)
    L <- loss(b) + lambda * t(b) %*% b

    maxIter <- 10
    tol <- 1e-6

    for (i in 1:maxIter) {
        w <- pmax(dinvlink(X %*% b), 1e-6)
        W <- Diagonal(x = as.numeric(w))
        H <- t(X) %*% W %*% X + 2 * lambda * diag(rep(1, d))
        grad <- gradloss(b)

        deltab <- newton_step(grad, H)

        step <- backtrack_linesearch(loss, gradloss, b, deltab, 0.3, 0.2)

        b <- b + step * deltab
        # b <- b - solve(H, grad)
        L_new <- loss(b) + lambda * t(b) %*% b

        if (as.logical(abs(L - L_new) < tol)) {
            # print(sprintf("No. iterations: %d", i))
            break
        }
        L <- L_new
    }
    return(b)
}


#defining an object class for a logistic regression model
logistic_model <- setRefClass("logistic_model",
                                 fields = c(
                                   X = "matrix",
                                   y = "numeric",
                                   coeffs = "numeric",
                                   prob = "numeric"))

#to initialise provide a design matrix and output label
#then uses the logistic regression implementation to find the coefficients
#uses the coefficients to find p(y=1) (probability of signal)

logistic_model$methods(
  initialize = function(X_train, y_train, lambda=1e-6){
    X_train <- as.matrix(X_train)
    y_train <- as.numeric(y_train)

    .self$X <- X_train
    .self$y <- y_train

    b <- logistic_reg(X_train, y_train, lambda=lambda)
    .self$coeffs <- as.numeric(b)
    .self$prob <- as.numeric(logistic(X_train %*% b))
  },

  #also defines a method to use this model to predict class of new points
  predict  = function(X_test){
    return(logistic(as.matrix(X_test) %*% .self$coeffs))
  }
)
