library(ggplot2)
source("utility_funcs.r")

#' Compute receiver operating characteristic (ROC) curve
#'
#' @param y response vector
#' @param p_hat probability outputs from model
#' @return list of false positive and true positive rates at different thresholds
compute_roc <- function(y, p_hat) {
    # FP = false positive rate = no. false positives / number of negatives
    # N <- 50
    # s <- 1/N
    # FP <- rep(NA, N)
    # TP <- rep(NA, N)
    # for (t in 1:N) {
    #   FP[t] <- sum(p_hat > t*s & y == 0) / sum(y == 0)
    #   TP[t] <- sum(p_hat > t*s & y == 1) / sum(y == 1)
    # }
    # return(list(FP, TP))

    #instead of thresholding at equally spaced thresholds, do it every n samples
    #because otherwise the points on the ROC are very uneven
    #order samples by p_hat
    y <- y[order(p_hat)]
    p_hat <- p_hat[order(p_hat)]

    #take the 1st, and last p-value, and 20 evenly spaced in the middle
    ints <- c(seq(1, length(p_hat), length(p_hat)/20), length(p_hat))
    p_ints <- p_hat[round(ints)]

    #for each p value, calculate the false positive rate and false negative rate
    FP <- rep(NA, length(p_ints))
    TP <- rep(NA, length(p_ints))
    for(i in 1:length(p_ints)){
        FP[i] <- sum(p_hat > p_ints[i] & y == 0) / sum(y == 0)
        TP[i] <- sum(p_hat > p_ints[i] & y == 1) / sum(y == 1)
    }
    return(list(FP, TP))
}

#' Calculate the area under the ROC curve (AUC) as a metric of performance
#'
#' @param FP false positive rate (vector)
#' @param TP true positive rate (vector)
#' @return AUC estimate (scalar)
compute_AUC <- function(FP, TP) {
    #using trapezoidal rule to find integral over our points
    AUC <- 0
    FP <- rev(FP)
    TP <- rev(TP)
    for(i in 1:(length(FP)-1)){
      h <- FP[i+1] - FP[i]
      AUC <- AUC + (h/2 * (TP[i] + TP[i+1]))
    }

    #f <- splinefun(FP, TP)
    #AUC <- integrate(f, 0, 1)
    return(AUC)
}

#' Plot ROC curve for particular model
#'
#' @param y response vector
#' @param p_hat model output probabilities
plot_roc <- function(y, p_hat, modelnum="k j") {
    # ROC plots FP vs TP. To get curve we vary the threshold
    # FP = false positive rate = no. false positives / number of negatives
    roc <- compute_roc(y, p_hat)
    FP <- roc[[1]]
    TP <- roc[[2]]
    auc <- compute_AUC(FP, TP)
    plot(0:1, 0:1, type="l", lty=2, col="red", xlab="False Positive Rate",
         ylab="True Positive Rate", main=paste("Model:", modelnum, ", AUC:", round(auc,2)))
    lines(FP, TP)
    legend("bottomright", legend=c("Chance", "Logistic Regression"), col=2:1, lty=2:1, cex=0.8)
}

#' plots of parameter values by fold to compare and check consistency of models
#' @param b matrix of coefficients (dxK)
plot_coefs <- function(b) {
    d <- nrow(b)
    K <- ncol(b)
    min_y <- min(b, na.rm=TRUE)
    max_y <- max(b, na.rm=TRUE)
    kcolours <- c("red", "blue", "green", "black", "orange")
    plot(c(1,d), c(min_y, max_y), type="n")
    for (k in 1:K) {
        points(1:d, b[,k], col=kcolours[k])
    }
}

ROC_curve <- setRefClass("ROC_curve",
                                 fields = c(
                                   thresholds = "numeric",
                                   FP = "numeric",
                                   TP = "numeric",
                                   auc = "numeric"))

ROC_curve$methods(
  initialize = function(y, p_hat){
    #we want just a set of thresholds that vary FP and FN rate from 0 to 1

    #order samples by p_hat
    y <- y[order(p_hat)]
    p_hat <- p_hat[order(p_hat)]

    #find 50 evenly spaced values from the list of probabilities, including 0, lowest, highest and 1
    indices <- 1 + 0:47 * (length(y)-1)/47
    .self$thresholds <- c(0, p_hat[round(indices)], 1)

    #for each threshold, calculate the false positive rate and false negative rate
    FPr <- rep(NA, 50)
    TPr <- rep(NA, 50)
    for(i in 1:50){
        FPr[i] <- sum(p_hat > thresholds[i] & y == 0) / sum(y == 0)
        TPr[i] <- sum(p_hat > thresholds[i] & y == 1) / sum(y == 1)
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

  #define a method to plot the roc (if add is true, only the new roc line is plotted)
  plot_curve = function(title="auc", add=FALSE){
    if(title=="auc"){
      title <- paste("AUC:", round(.self$auc,2))
    }
    if(add==FALSE){
      plot(0:1, 0:1, type="l", lty=2, col="red", xlab="False Positive Rate",
           ylab="True Positive Rate", main=title)
      legend("bottomright", legend=c("Chance", "Logistic Regression"), col=2:1, lty=2:1)
    }
    lines(.self$FP, .self$TP)
  }
)


AMS_data <- setRefClass("AMS_data",
                        fields = c(
                          y = "numeric",
                          prob = "numeric",
                          weights = "numeric",
                          sum_w = "ANY",
                          thresholds = "numeric",
                          ams = "numeric",
                          max_ams = "numeric",
                          max_thresh = "numeric"))
AMS_data$methods(
  #to initialise provide a list of true labels, probabilities P(y="s"), and scaled weights
  initialize = function(y, prob, weights, sum_w=NULL){
    .self$y <- as.numeric(y)
    .self$prob <- as.numeric(prob)
    .self$weights <- as.numeric(weights)
    .self$sum_w <- if (!is.null(sum_w)) as.numeric(sum_w)
    .self$thresholds <- seq(0, 1, length.out = 30)
  },

  calc_ams = function(){
    #for each threshold, calculate b, s, and the AMS
    n <- length(thresholds)
    AMS <- rep(NA, n)

    for(i in 1:n){
      #use the threshold to make a classification
      y_pred <- prob >= thresholds[i]
      AMS[i] <- calculate_ams_partition(y, y_pred, weights, sum_w)
    }
    .self$ams <- AMS
  },

  get_max_thresh = function(){
    .self$max_thresh <- thresholds[which.max(ams)]
  },

  plot_ams = function(lgd=NULL, add=FALSE){
    #find the maximum ams threshold
    .self$calc_ams()
    .self$max_ams <- max(ams)
    .self$get_max_thresh()

    #plot the ams at different thresholds, with a line at the best
    if(add==FALSE){
      plot(thresholds, ams, type="l", col="steelblue", main="AMS at different decision thresholds")
      abline(v=max_thresh, lty=2)
      if (is.null(lgd)) {lgd <- paste0("Max AMS at p=", round(max_thresh, 2))}
      legend("bottomleft", legend=lgd, lty=2)
    }
    lines(thresholds, ams)
  }
)
