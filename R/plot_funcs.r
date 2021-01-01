library(ggplot2)
#source("utility_funcs.r") 

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

# want to plot a curve for each fold
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

# plots of parameter values by fold to compare.
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

# plot scatter matrix
# library(psych)
# pairs.panels(head(X[no_missing, 1:10], 1000), 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = TRUE # show correlation ellipses
#              )