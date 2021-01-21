#' @importFrom CVXR Variable Problem
fit_l1_logreg <- function(X, y, C=1) {
  nb <- ncol(X)
  b <- Variable(nb, 1)
  obj <- -sum(X[y == 0, ] %*% b) - sum(logistic(-X %*% b))
  constr <- list(abs(b[2:nb]) <= C)
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob)
  b_res_con <- result$getValue(b)
}

logisticf <- function(x) {
  1/(1 + exp(-x))
}

#defining an object class for a logistic regression model
#' @import methods
logistic_l1_model <- setRefClass("logistic_l1_model",
                              fields = c(
                                X = "matrix",
                                y = "numeric",
                                coeffs = "numeric",
                                prob = "numeric"))

#to initialise provide a design matrix and output label
#then uses the logistic regression implementation to find the coefficients
#uses the coefficients to find p(y=1) (probability of signal)
#' @import methods
logistic_l1_model$methods(
  initialize = function(X_train, y_train, C=1){
    X_train <- as.matrix(X_train)
    y_train <- as.numeric(y_train)

    .self$X <- X_train
    .self$y <- y_train

    b <- fit_l1_logreg(X_train, y_train, C=C)
    .self$coeffs <- as.numeric(b)
    .self$prob <- as.numeric(logisticf(X_train %*% b))
  },

  #also defines a method to use this model to predict class of new points
  predict = function(X_test){
    return(logisticf(as.matrix(X_test) %*% .self$coeffs))
  }
)
