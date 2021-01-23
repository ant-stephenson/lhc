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
#' @import RColorBrewer
#' @export
plot_rocs <- function(rocs, title="ROC curves", ...){
  n <- length(rocs)
  colours <- colorRampPalette(brewer.pal(8, 'Blues')[2:8])(n) #gets some nice blues
  plot(0:1, 0:1, type="l", lty=2, xlab="False Positive Rate",
       ylab="True Positive Rate", main=title, ...)
  legend("bottomright", legend=c("Chance", "Logistic Regression"),
         col=c("black", colours[n]), lty=2:1)
  for(i in 1:n){
    lines(rocs[[i]]$FP, rocs[[i]]$TP, col=colours[i])
  }
}

#' define a function to plot multiple ams objects on the same axes
#' @param amss list of ams objects
#' @param title str title to give plot
#' @import RColorBrewer
#' @export
plot_amss <- function(amss, title="AMS data", min.max=TRUE, ...){
  lapply(amss, function(x) x$calc_ams())
  y_max <- max(sapply(amss, function(x) max(x$ams)))
  n <- length(amss)
  colours <- colorRampPalette(brewer.pal(8, 'Blues')[2:8])(n) #gets some nice blues
  ams <- amss[[1]]

  # calculate the minimum bounding curve of the curves in amss
  if (min.max) {
    ams_vec <- sapply(amss, function(x) x$ams)
    min_ams_obj <- ams
    min_ams <- pmax(apply(ams_vec, 1, min), 1e-1)
    min_ams_obj$ams <- min_ams
    min_max_thresh <- min_ams_obj$get_max_thresh()

    title <- sprintf("AMS plot with max-min (over folds) threshold at t=%.3f", min_max_thresh)
  }

  plot(ams$thresholds, ams$ams, type="l", col=colours[1], main=title,
       xlab="Decision Threshold", ylab="AMS", ylim=c(0, y_max), ...)

  for(i in 2:n){
    ams <- amss[[i]]
    lines(ams$thresholds, ams$ams, type="l", col=colours[i])
    if (!min.max) {
      abline(v=ams$max_thresh, lty=2)
    }
  }

  if (min.max) {
    abline(v=min_max_thresh, lty=2)
  }
  legend=legend("bottomleft", legend=1:n, fill=colours)
}

#' function to wrap figure saving
#' @param plot_func partially called function of no arguments to generate plot
#' @param filepath string file path to save to
#' @param filetype type to save as (function name)
#' @export
save_fig <- function(plot_func, filepath, filetype=pdf) {
  filetype(file=filepath)
  plot_func()
  dev.off() -> .
}

#' Density plots of variables
#'
#' Given a matrix plot the density of the listed variables using ggplot2 facet_wrap.
#'
#' If the matrix X is more than 10,000 samples, a random 10,000 samples will be selected
#' to keep the amount of data plotted reasonable.
#'
#' @param X An nxd matrix with samples as rows and features as columns.
#' @param variables Optional vector of column names to be plotted
#' @param lables Optional vector of class labels to view distributions by class
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @export
plot_distributions <- function(X, variables=NULL, labels=NULL){

  #if labels not provided just use 1s
  if(is.null(labels)){
    labels <- rep(1, nrow(X))
  }

  #if X is large, take a sample
  if(nrow(X) > 10000){
    idx <- sample(1:nrow(X), 10000)
    X <- X[idx,]
    labels <- labels[idx]
  }

  #if columns specified, select just these
  if(!is.null(variables)){
    X <- X[, variables]
  }

  #unpivot (make a long format)
  plot_data <- cbind(rownames(X), labels, as.data.frame(X)) %>%
    pivot_longer(cols=c(-1, -2), values_to = "Value", names_to="Variable")

  #create plot
  ggplot(plot_data, aes(Value, colour=labels)) +
    geom_density() +
    theme_minimal() +
    theme(legend.position = "none") +
    facet_wrap(vars(Variable), scales="free")
}
