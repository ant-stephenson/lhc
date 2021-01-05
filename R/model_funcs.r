source("utility_funcs.R")

#' Fit a logistic regression model by IRWLS - same thing glm does
#'
#' @param X covariate matrix
#' @param y response vector
#' @param r [Optional] weight vector
#' @param lambda [Optional] L2 regularisation parameter
#' @return b vector of coefficients
library(Matrix)
logistic_reg <- function(X, y, r=NULL, lambda = 0) {
    invlink <- logistic
    dinvlink <- function(x)  exp(x)/(1 + exp(x))^2
    loss <- function(y, eta) log(1 + exp((-1)^y*eta))

    d <- ncol(X)
    n <- nrow(X)

    if (is.null(r)) {r <- rep(1,n)}

    b <- matrix(0, d, 1)
    eta <- matrix(0, n, 1)
    mu <- invlink(eta)
    L <- sum(loss(y, eta)) + lambda * t(b) %*% b

    maxIter <- 10
    tol <- 1e-6

    for (i in 1:maxIter) {
        w <- pmax(dinvlink(eta), 1e-6)
        W <- Diagonal(x = as.numeric(w) * r)
        H <- t(X) %*% W %*% X + 2 * lambda * diag(rep(1, d))
        grad <- -t(X) %*% (y - mu) + 2 * lambda * b

        #an option that does not compute H^-1
        #marginally quicker..
        H_svd <- svd(H)
        H_inv <- H_svd$v %*% diag(1/H_svd$d) %*% t(H_svd$u)
        b <- b - H_inv %*% grad
        # b <- b - solve(H, grad)
        eta <- X %*% b
        mu <- invlink(eta)
        L_new <- sum(loss(y, eta)) + lambda * t(b) %*% b

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
  predict = function(X_test){
    return(logistic(as.matrix(X_test) %*% .self$coeffs))
  }
)
