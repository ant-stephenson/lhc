#' Generate distinct colours
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


#' Calculate average AUC
#' @param rocs a list of ROC_curve objects
#' @return average AUC
average_auc <- function(rocs){
  lapply(rocs, function(x) x$calc_auc())
  auc <- mean(sapply(rocs, function(x) x$auc))
  return(auc)
}


#' Plot Multiple ROC curves
#' @param rocs a list of ROC_curve objects
#' @param title str title to give plot
#' @param info additional info to add to the default title
#' @import RColorBrewer
#' @export
plot_rocs <- function(rocs, title=NULL, info="", ...){
  n <- length(rocs)
  auc <- average_auc(rocs)
  colours <- colorRampPalette(brewer.pal(8, 'Blues')[2:8])(n) #gets some nice blues

  if(is.null(title)){
    title <- sprintf("%s ROC curves with mean AUC of %.3f", info, auc)
  }

  plot(0:1, 0:1, type="l", lty=2, xlab="False Positive Rate",
       ylab="True Positive Rate", main=title, ...)

  for(i in 1:n){
    lines(rocs[[i]]$FP, rocs[[i]]$TP, col=colours[i])
  }

  legend("bottomright", legend=c("Chance", paste0("Models 1-", n)),
         col=c("black", colours[n]), lty=2:1, cex=0.8)
}

#' Calculate optimal threshold
#' @param amss a list of AMS_data objects
#' @return optimal decision threshold considering all curves#
ams_threshold <- function(amss){
  # calculate the minimum bounding curve of the curves in amss
  ams_vec <- sapply(amss, function(x) x$ams)
  min_ams_obj <- amss[[1]]
  min_ams <- pmax(apply(ams_vec, 1, min), 1e-1)
  min_ams_obj$ams <- min_ams
  min_max_thresh <- min_ams_obj$get_max_thresh()
  return(min_max_thresh)
}

#' define a function to plot multiple ams objects on the same axes
#' @param amss list of ams objects
#' @param title str title to give plot, if null a default title is generated
#' @param info additional info to add to the default title
#' @import RColorBrewer
#' @export

plot_amss <- function(amss, title=NULL, info="", min.max=TRUE, ...){
  lapply(amss, function(x) x$calc_ams())
  y_max <- max(sapply(amss, function(x) max(x$ams)))
  n <- length(amss)
  colours <- colorRampPalette(brewer.pal(8, 'Blues')[2:8])(n) #gets some nice blues
  ams <- amss[[1]]

  min_max_thresh <- ams_threshold(amss)
  if(is.null(title)){
    title <- sprintf("%s AMS plots with threshold at %.3f", info, min_max_thresh)
  }

  plot(ams$thresholds, ams$ams, type="l", col=colours[1], main=title,
       xlab="Decision Threshold", ylab="AMS", ylim=c(0, y_max), ...)

  for(i in 2:n){
    ams <- amss[[i]]
    lines(ams$thresholds, ams$ams, type="l", col=colours[i])
  }
  abline(v=min_max_thresh, lty=2)

  legend("bottomleft", legend=c("Optimal threshold", paste0("Models 1-", n)),
                col=c("black", colours[n]), lty=2:1, cex=0.8)
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
