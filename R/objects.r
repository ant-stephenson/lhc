#' Logistic model object class
#'
#' A reference class object to hold a logistic regression model.
#' Can be used to initialize and train a model using the logistic_reg function,
#' then use the model to predict the classes of new samples.
#' Sample classes must be 0 or 1.
#' Predict returns the probability that each new sample is class 1,
#' and a decision threshold should be subsequently applied.
#'
#' @name logistic_model
#' @import methods
#' @export logistic_model
#' @exportClass logistic_model
#'
#' @field X An nxd matrix with samples as rows and features as columns
#' @field y A length-n vector of 0s and 1s indicating true sample classes
#' @field coeffs A length-d vector of model coefficients
logistic_model <- setRefClass("logistic_model",
                              fields = c(
                                X = "matrix",
                                y = "numeric",
                                coeffs = "numeric"
                              ),
                              methods = list(
                                initialize = function(X, y){
                                  "Provide X and y and the coeffs field will be calculated using logistic_reg"
                                  .self$X <- as.matrix(X)
                                  .self$y <- as.numeric(y)
                                  .self$coeffs <- as.numeric(logistic_reg(as.matrix(X), as.numeric(y)))
                                },
                                predict = function(X_test){
                                  "Provide a matrix of new samples and a vector of P(y=1) is returned"
                                  prob <- logisticf(as.matrix(X_test) %*% coeffs)
                                  return(as.numeric(prob))
                                }
                              )
)

ROC_data <- setRefClass("ROC_data",
                        fields = c(
                          thresholds = "numeric",
                          FP = "numeric",
                          TP = "numeric",
                          auc = "numeric"))
ROC_data$methods(
  initialize = function(y, prob){
    #we want to find a set of thresholds that vary FP and FN rate from 0 to 1
    #find n evenly spaced values from the list of probabilities, including the 0 and 1
    prob_ordered <- prob[order(prob)]
    n <- 30
    indices <- 1 + 0:(n-3) * (length(y)-1)/(n-3)
    thresh <- c(0, prob_ordered[round(indices)], 1)

    #for each threshold, calculate the FR, FN, s and b
    FPr <- rep(NA, n)
    TPr <- rep(NA, n)
    for(i in 1:n){
      FPr[i] <- sum(prob > thresh[i] & y == 0) / sum(y == 0)
      TPr[i] <- sum(prob > thresh[i] & y == 1) / sum(y == 1)
    }

    .self$FP <- FPr
    .self$TP <- TPr
    .self$auc <- 0
  },

  calc_auc = function(){
    if(.self$auc == F){
      #using the trapezoidal rule to integrate the roc curve to find the auc
      #since we only have points not a curve I think this is approach is sufficent
      #note that the points go from right to left on the curve (so FP[i] > FP[i+1])
      AUC <- 0
      for(i in 1:(length(FP)-1)){
        h <- FP[i] - FP[i+1]
        AUC <- AUC + (h/2 * (TP[i] + TP[i+1]))
      }
      .self$auc <- AUC
      return(AUC)
    }
  },

  #define a method to plot the roc
  plot_curve = function(){
    .self$calc_auc()
    plot(0:1, 0:1, type="l", lty=2, xlab="False Positive Rate",
         ylab="True Positive Rate", main=paste("AUC:", round(.self$auc,2)))
    lines(.self$FP, .self$TP, col="steelblue")
    legend("bottomright", legend=c("Model", "Chance"), col=c("steelblue", 1), lty=1:2)
  }
)


AMS_data <- setRefClass("AMS_data",
                        fields = c(
                          y = "numeric",
                          prob = "numeric",
                          weights = "numeric",
                          thresholds = "numeric",
                          ams = "numeric"))
AMS_data$methods(
  #to initialise provide a list of true labels, probabilities P(y="s"), and scaled weights
  initialize = function(y, prob, weights){
    .self$y <- as.numeric(y)
    .self$prob <- as.numeric(prob)
    .self$weights <- as.numeric(weights)
    .self$thresholds <- seq(0, 1, length.out = 30)
  },

  calc_ams = function(){
    #for each threshold, calculate b, s, and the AMS
    n <- length(thresholds)
    AMS <- rep(NA, n)

    for(i in 1:n){
      #use the threshold to make a classification
      y_pred <- prob >= thresholds[i]
      AMS[i] <- ams_metric(y, y_pred, weights)
    }
    .self$ams <- AMS
  },

  plot_ams = function(){
    #find the maximum ams threshold
    .self$calc_ams()
    max_ams <- max(ams)
    max_thresh <- thresholds[which.max(ams)]

    #plot the ams at different thresholds, with a line at the best
    plot(thresholds, ams, type="l", col="steelblue", main="AMS at different decision thresholds")
    abline(v=max_thresh, lty=2)
    legend("topright", legend=paste0("Max AMS at p=", round(max_thresh, 2)), lty=2)
  }
)
