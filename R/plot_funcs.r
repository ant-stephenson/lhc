#' Function to generate colours that are quite distinct
#' @param ncolours number of colours we want
#' @return vector of hex colours
#' @export
generate_colours <- function(ncolours) {
  colours <- vector("list", ncolours)
  for (i in 1:ncolours) {
    h <- i/ncolours
    s <- 0.8 + 0.1 * runif(1)
    v <- 0.6 + 0.1 * runif(1)
    colours[i] <- hsv(h,s,v)
  }
  return(unlist(colours))
}

#If we're going for OOP approach for ROC, use ROC_curve$new(y_test, prob)
#' Compute receiver operating characteristic (ROC) curve
#'
#' @param y response vector
#' @param p_hat probability outputs from model
#' @return list of false positive and true positive rates at different thresholds
# compute_roc <- function(y, p_hat) {
#     #threshold every n samples
#     #because otherwise the points on the ROC are very uneven
#     #order samples by p_hat
#     y <- y[order(p_hat)]
#     p_hat <- p_hat[order(p_hat)]
#
#     #take the 1st, and last p-value, and 20 evenly spaced in the middle
#     ints <- c(seq(1, length(p_hat), length(p_hat)/20), length(p_hat))
#     p_ints <- p_hat[round(ints)]
#
#     #for each p value, calculate the false positive rate and false negative rate
#     FP <- rep(NA, length(p_ints))
#     TP <- rep(NA, length(p_ints))
#     for(i in 1:length(p_ints)){
#         FP[i] <- sum(p_hat > p_ints[i] & y == 0) / sum(y == 0)
#         TP[i] <- sum(p_hat > p_ints[i] & y == 1) / sum(y == 1)
#     }
#     return(list(FP, TP))
# }

#If we're going for OOP approach for ROC, use ROC_curve$new(y_test, prob); roc$calc_auc()
#' Calculate the area under the ROC curve (AUC) as a metric of performance
#'
#' @param FP false positive rate (vector)
#' @param TP true positive rate (vector)
#' @return AUC estimate (scalar)
# compute_AUC <- function(FP, TP) {
#     #using trapezoidal rule to find integral over our points
#     AUC <- 0
#     FP <- rev(FP)
#     TP <- rev(TP)
#     for(i in 1:(length(FP)-1)){
#       h <- FP[i+1] - FP[i]
#       AUC <- AUC + (h/2 * (TP[i] + TP[i+1]))
#     }
#
#     #f <- splinefun(FP, TP)
#     #AUC <- integrate(f, 0, 1)
#     return(AUC)
# }

#If we're going for OOP approach for ROC, use ROC_curve$new(y_test, prob); roc$plot_curve()
#' Plot ROC curve for particular model
#'
#' @param y response vector
#' @param p_hat model output probabilities
# plot_roc <- function(y, p_hat, modelnum="k j") {
#     # ROC plots FP vs TP. To get curve we vary the threshold
#     # FP = false positive rate = no. false positives / number of negatives
#     roc <- compute_roc(y, p_hat)
#     FP <- roc[[1]]
#     TP <- roc[[2]]
#     auc <- compute_AUC(FP, TP)
#     plot(0:1, 0:1, type="l", lty=2, col="red", xlab="False Positive Rate",
#          ylab="True Positive Rate", main=paste("Model:", modelnum, ", AUC:", round(auc,2)))
#     lines(FP, TP)
#     legend("bottomright", legend=c("Chance", "Logistic Regression"), col=2:1, lty=2:1, cex=0.8)
# }

#Is this used?
#' plots of parameter values by fold to compare and check consistency of models
#' @param b matrix of coefficients (dxK)
# plot_coefs <- function(b) {
#     d <- nrow(b)
#     K <- ncol(b)
#     min_y <- min(b, na.rm=TRUE)
#     max_y <- max(b, na.rm=TRUE)
#     kcolours <- c("red", "blue", "green", "black", "orange")
#     plot(c(1,d), c(min_y, max_y), type="n")
#     for (k in 1:K) {
#         points(1:d, b[,k], col=kcolours[k])
#     }
# }


#' define a function to plot multiple roc objects on the same axes
#' @param rocs list of roc objects
#' @param title str title to give plot
#' @export
plot_rocs <- function(rocs, title="ROC curves"){
  n <- length(rocs)
  colours <- brewer.pal(n+1, 'Blues')[2:(n+1)] #gets some nice blues but skips the lightest one
  plot(0:1, 0:1, type="l", lty=2, xlab="False Positive Rate",
       ylab="True Positive Rate", main=title)
  legend("bottomright", legend=c("Chance", "Logistic Regression"),
         col=c("black", colours[n]), lty=2:1)
  for(i in 1:n){
    lines(rocs[[i]]$FP, rocs[[i]]$TP, col=colours[i])
  }
}

#' define a function to plot multiple ams objects on the same axes
#' @param amss list of ams objects
#' @param title str title to give plot
#' @export
plot_amss <- function(amss, title="AMS data"){
  lapply(amss, function(x) x$calc_ams())
  y_max <- max(sapply(amss, function(x) max(x$ams)))
  n <- length(amss)
  colours <- brewer.pal(n+1, 'Blues')[2:(n+1)] #gets some nice blues but skips the lightest one
  ams <- amss[[1]]
  plot(ams$thresholds, ams$ams, type="l", col=colours[1], main=title,
       xlab="Decision Threshold", ylab="AMS", ylim=c(0, y_max))

  for(i in 2:n){
    ams <- amss[[i]]
    lines(ams$thresholds, ams$ams, type="l", col=colours[i])
  }
}
