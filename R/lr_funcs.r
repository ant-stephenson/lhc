#library(purrr)
#source("utility_funcs.R")

#' Compute newton step
#' @param grad gradient vector
#' @param H Hessian matrix
#' @return step to take
#' @import purrr
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
#'
#' @import purrr
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

#' x contains lambda (dual var) as well
#' @param f objective function
#' @param dualf dual objective function
#' @param gradf residual vector
#' @param Hf Hessian for residual
#' @param x primal-dual point
#' @param m number of inequality constraints
#' @param mu interior-point step parameter
#' @param eps tolerance for problem
#' @param eps_feas tolerance for feasibility of primal-dual points
interior_point_fit <- function(f, dualf, gradf, Hf, x, m, mu=10, eps=1e-6, eps_feas=1e-6) {
  conv_cond <- function(x, tt, eta) {
    all(abs(gradf(x, tt)) <= eps_feas) & (eta <= eps)
  }

  eta <- f(x) - dualf(x)
  tt <- m/eta

  max_iter <- 100
  iter <- 0
  while (!conv_cond(x, tt, eta) & iter < max_iter) {
    print(eta)
    iter <- iter + 1
    tt <- mu * m / eta
    delta <- newton_step(gradf(x, tt), Hf(x))
    step <- backtrack_linesearch(f, partial(gradf, tt=tt) , x, delta, 0.3, 0.2)
    x <- x + step * delta
    eta <- f(x) - dualf(x)
  }
  print(iter)
  return(x)
}

#' @import purrr
l1_logistic_reg <- function(X, y, C) {
  invlink <- logisticf
  dinvlink <- function(x) exp(-x)/(1+exp(-x))^2
  loss <- function(b) {
    nx <- length(b)/3
    sum(log(1 + exp((-1)^y*(X %*% b[1:nx]))))
  }
  grad <- function(b, tt) {
    nx <- length(b)/3
    m <- 2 * nx
    lambda <- b[-(1:nx)]
    g <- -t(X) %*% (y - invlink(X %*% b[1:nx])) + sum(lambda)
    g <- rbind(g, -t(t(lambda * c(-b[1:nx]-C, b[1:nx]-C) - 1/tt)))
  }
  hess <- function(b) {
    nx <- length(b)/3
    m <- 2 * nx
    lambda <- b[-(1:nx)]
    w <- pmax(dinvlink(X %*% b[1:nx]), 1e-6)
    W <- Diagonal(x = as.numeric(w))
    H <- t(X) %*% W %*% X
    ones <- rbind(diag(rep(1,nx)), diag(rep(1,nx)))
    H <- cbind(H, t(ones))
    extra_rows <- cbind(-ones * lambda, -diag(c(-b[1:nx] - C, b[1:nx] - C)))
    H <- rbind(H, extra_rows)
  }
  dual <- function(b) {
    nx <- length(b)/3
    m <- 2 * nx
    lambda <- b[-(1:nx)]
    return(loss(b) + sum(lambda * (abs(b[1:nx]) - C)))
  }
  d <- ncol(X)
  n <- nrow(X)

  b <- matrix(c(rep(0, d), rep(1, 2*d)), 3*d, 1)

  b <- interior_point_fit(loss, dual, grad, hess, b, 2*d, mu=10, eps=1e-6, eps_feas=1e-6)
}

#' Fit a logistic regression model by IRWLS - same thing glm does
#'
#' @param X covariate matrix
#' @param y response vector
#' @param r [Optional] weight vector
#' @param lambda [Optional] L2 regularisation parameter
#' @return b vector of coefficients
#' @import Matrix
logistic_reg <- function(X, y, lambda = 0) {
    invlink <- logisticf
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
#' @import methods
logistic_model <- setRefClass("logistic_model",
                                 fields = c(
                                   X = "matrix",
                                   y = "numeric",
                                   coeffs = "numeric",
                                   prob = "numeric"))

#to initialise provide a design matrix and output label
#then uses the logistic regression implementation to find the coefficients
#uses the coefficients to find p(y=1) (probability of signal)
#' @import methods
logistic_model$methods(
  initialize = function(X_train, y_train, lambda=1e-6){
    X_train <- as.matrix(X_train)
    y_train <- as.numeric(y_train)

    .self$X <- X_train
    .self$y <- y_train

    b <- logistic_reg(X_train, y_train, lambda=lambda)
    .self$coeffs <- as.numeric(b)
    .self$prob <- as.numeric(logisticf(X_train %*% b))
  },

  #also defines a method to use this model to predict class of new points
  predict  = function(X_test){
    return(logisticf(as.matrix(X_test) %*% .self$coeffs))
  }
)

# X <- model.matrix(~1 + x1 + x2, data.frame(x1=rnorm(1000), x2=rnorm(1000)))
# b <- abs(rnorm(3)) + 1
# y <- rbinom(1000, 1, logisticf(X %*% b))
# mod <- l1_logistic_reg(X, y, 1)
# mod2 <- logistic_reg(X, y)
