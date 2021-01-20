#' Logistic model object class
#'
#' This reference class object fits a binary classification model,
#' using the logistic_reg function. The model can be used to
#' predict the classes of new samples. The sample classes must be 0 or 1,
#' and the prediction returns the estimated probabilities that each sample is class 1.
#' A decision threshold should be subsequently.
#'
#' @name logistic_model
#' @import methods
#' @export logistic_model
#' @exportClass logistic_model
#'
#' @field X An nxd matrix with samples as rows and features as columns.
#' @field y A length-n vector of 0s and 1s indicating true sample classes.
#' @field coeffs A length-d vector of model coefficients.
#' @field lambda regularization parameter
logistic_model <- setRefClass("logistic_model",
  fields = c(
    X = "matrix",
    y = "numeric",
    coeffs = "numeric",
    lambda = "numeric"
  ),
  methods = list(
    initialize = function(X, y, lambda=1e-6){
      "Provide X and y and the coeffs field will be calculated using logistic_reg"
      .self$X <- as.matrix(X)
      .self$y <- as.numeric(y)
      .self$coeffs <- as.numeric(logistic_reg(as.matrix(X), as.numeric(y), lambda=lambda))
    },
    predict = function(X_test){
      "Provide a matrix of new samples and a vector of P(y=1) is returned"
      prob <- logisticf(as.matrix(X_test) %*% coeffs)
      return(as.numeric(prob))
    }
  )
)

#' ROC curve object class
#'
#' This reference class object is used to plot a Receiver Operating Characteristic curve.
#' An ROC curve is a performance measure of a classification model,
#' created by plotting the true positive rate (TPR) against the false positive rate (FPR)
#' as the decision threshold is varied. The object finds suitable thresholds,
#' calculates FPR and TPR at each, and can calculate the Area Under the Curve (AUC).
#' A vector of true sample classifications (0 or 1) and a vector of estimated probabilities
#' from a model are needed to initialise.
#'
#' @name ROC_curve
#' @import methods
#' @export ROC_curve
#' @exportClass ROC_curve
#'
#' @field thresholds A vector of 30 decision thresholds.
#' @field FP A vector of the false positive rate at each threshold.
#' @field TP A vector of the true positive rate at each threshold.
#' @field auc A numeric that is the area under the ROC curve.
ROC_curve <- setRefClass("ROC_curve",
  fields = c(
    thresholds = "numeric",
    FP = "numeric",
    TP = "numeric",
    auc = "numeric"
  ),
  methods = list(
    initialize = function(y, prob){
      "Provide sample labels and probabilites, and the FPR and TPR are calculated at 30 decision thresholds"
      #find n thresholds that vary FP and FN rate from 0 to 1
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

      #add results to the object fields
      .self$thresholds <- thresh
      .self$FP <- FPr
      .self$TP <- TPr
      .self$auc <- 0
    },
    calc_auc = function(){
      "If the AUC has not already be calculated, this calls the calculation."
      if(.self$auc == F){
        #use the trapezoidal rule to integrate the ROC curve to find the AUC
        #since we have n points not a curve this is sufficient
        AUC <- 0
        for(i in 1:(length(FP)-1)){
          h <- FP[i] - FP[i+1]
          AUC <- AUC + (h/2 * (TP[i] + TP[i+1]))
        }
        .self$auc <- AUC
        return(AUC)
      }
    },
    plot_curve = function(){
      "Plot the ROC curve."
      .self$calc_auc()
      plot(0:1, 0:1, type="l", lty=2, xlab="False Positive Rate",
           ylab="True Positive Rate", main=paste("AUC:", round(.self$auc,2)))
      lines(.self$FP, .self$TP, col="steelblue")
      legend("bottomright", legend=c("Model", "Chance"), col=c("steelblue", 1), lty=1:2)
    }
  )
)

#' AMS data object class
#'
#' This reference class object is used to store the AMS metric of a classification
#' model at different decision thresholds. AMS is a performance measure which
#' includes the sample weightings and is defined by the Higgs Boson Kaggle Competition.
#' A vector of true sample classifications (0 or 1), a vector of estimated probabilities
#' from a model, and a vector of scaled sample weights are needed to initialise.
#'
#' @name AMS_data
#' @import methods
#' @export AMS_data
#' @exportClass AMS_data
#'
#' @field y A vector of true sample classifications (0 or 1),
#' @field prob A vector of the samples estimated probabilities from a model
#' @field weights A vector of scaled sample weights.
#' @field thresholds A vector of 30 decision thresholds.
#' @field ams A vector of the AMS metric at each threshold.
AMS_data <- setRefClass("AMS_data",
  fields = c(
    y = "numeric",
    prob = "numeric",
    weights = "numeric",
    thresholds = "numeric",
    ams = "numeric"
  ),
  methods = list(
    initialize = function(y, prob, weights){
      "Provide true sample lables, estimated probabilities, and sample weights. A vector of descision thresholds is initalised."
      .self$y <- as.numeric(y)
      .self$prob <- as.numeric(prob)
      .self$weights <- as.numeric(weights)
      .self$thresholds <- seq(0, 1, length.out = 30)
    },
    calc_ams = function(){
      "Calculate the AMS at each thresholds."
      n <- length(thresholds)
      AMS <- rep(NA, n)
      for(i in 1:n){
        #use each threshold to make a classification
        y_pred <- prob >= thresholds[i]
        #with these classifications calculate AMS
        AMS[i] <- calculate_ams_partition(y, y_pred, weights)
      }
      .self$ams <- AMS
    },
    plot_ams = function(){
      "Plot AMS against threshold."
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
)
