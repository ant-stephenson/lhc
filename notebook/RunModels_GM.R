
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(fs)
library(purrr)
library(lhc)

## ----load data---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# divide data into training kaggle set, and retain hold-out (before further cross-validation partitioning)
#source_path <- path_join(c(path_dir(getwd())))
#filename <- "atlas-higgs-challenge-2014-v2.csv"
#filepath <- path_join(c(source_path, "LHC_dump", "R", filename))
filepath <- list.files(path="~", pattern="atlas-higgs-challenge-2014-v2.csv",
                       full.names=T, recursive=T)

data <- import_data(filepath)
train_idx <- get_subset_idx(data$kaggle_s, "t")
y <- data$y[train_idx]
# need to use kaggle weignts to make AMS correct?!
kaggle_w <- data$kaggle_w[train_idx]
w <- data$kaggle_w[train_idx]
nj <- data$nj[train_idx]
G=3

# function just to loop over fitting procedure with various input parameters to faciliate running experiments over parameter space (e.g. a gridsearch)
run_models <- function(K, n_rbf, lambda, C, poly_order) {

  # modify features
  X <- data$X[train_idx, ]
  X <- reduce_features(X)
  X <- invert_angle_sign(X)

  #add polynomial transformations of order poly_order to our covariates
  if (poly_order > 1) {
    X <- poly_transform(X, poly_order)
  }

  ## ----init params-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  s <- avg_median_pairwise_distance(X)

  # K-Fold CV partitioning
  kI <- partition_data(length(y), K, random = TRUE)

  # number of rows and columns of data matrix
  n <- nrow(X)
  d <- ncol(X) + n_rbf

  ## ----jet/missing-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # get missing rows. separate by number of jets and presence of Higgs mass and fit separate models
  # find columns with features with any missing values for each number of jets: 0, 1, 2+ in combination with the presence (or absence) of the Higgs mass, defined by j=1,2,3, 4, 5, 6. i.e. j=1 => nj=0 & mH != -999, j=2 =>
  jet_cats <- 1:3
  features_to_rm <- set_features_to_rm(X, G, kI, nj)

  ## ----fitting-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # loop over folds
  #create lists to hold the k models and k roc curves
  models <- vector("list", G*K)
  rocs <- vector("list", G*K)
  ams_obj <- vector("list", G*K)

  ams <- rep(NA, length=G*K)
  auc <- rep(NA, length=G*K)

  for (mj in 1:G) {
    # loop over sets of jet number {0, 1, 2+} and mH presence/absence
    for (k in 1:K) {
      j <- jet_cats[mj]
      fit_row_idx <- kI != k & idx_jet_cat(nj, j)
      test_row_idx <- kI == k & idx_jet_cat(nj, j)

      # add r RBF centroid features, using the same reference centroids in training and testing sets
      if (n_rbf > 0) {
        rbf_centroids <- get_rbf_centroids(X[fit_row_idx, ], n_rbf)
        Xi <- rbf_centroids$"xi"
        Xtrain <- add_rbf_features(X[fit_row_idx, ], s, n_rbf, Xi=Xi)
        Xtest <- add_rbf_features(X[test_row_idx, ], s, n_rbf, Xi=Xi)
      } else {
        Xtrain <- X[fit_row_idx, ]
        Xtest <- X[test_row_idx, ]
      }

      # use group/fold to get valid columns
      col_idx <- get_valid_cols(colnames(X), features_to_rm, mj)

      # define the valid columns for this fold/group
      Xtrain <- as.matrix(Xtrain[, col_idx])
      Xtest <- as.matrix(Xtest[, col_idx])

      # standardize our datasets. Use training information to standardize the test set
      Xtest <- scale_dat(Xtest, Xtrain)
      Xtrain <- scale_dat(Xtrain, Xtrain)

      # get index for this model (defined order based on groups/folds)
      model_idx <- get_model_idx(mj, k, K)

      #fit a logistic regression model to the CV training data
      models[[model_idx]] <- logistic_model$new(X=Xtrain, y=y[fit_row_idx], lambda)

      #use it to predict the classifications of the test data
      p_hat <- models[[model_idx]]$predict(Xtest)

      #create an ROC curve object
      rocs[[model_idx]] <- ROC_curve$new(y[test_row_idx], p_hat)

      # create an AMS object
      w_this_partition <- w[test_row_idx]
      ams_obj[[model_idx]] <- AMS_data$new(y[test_row_idx], p_hat, w_this_partition)
    }
  }

  #get auc and ams from objects
  auc <- sapply(rocs, function(x) x$calc_auc())
  ams <- sapply(ams_obj, function(x) x$calc_ams())

  #get results for each group
  mean_auc <- mean_ams <- mad_ams <- mad_auc <- rep(NA, G)
  for (g in 1:G) {
    model_set <- ((g-1)*K+1):(g*K)
    mean_ams[g] <- mean(ams[,model_set])
    mean_auc[g] <- mean(auc[model_set])
    mad_ams[g] <- mean(apply(ams[, model_set], 1, mad))
    mad_auc[g] <- mad(auc[model_set])
  }

  result <- data.frame("K"=K, "G"=1:3, "n_rbf"=n_rbf, "lambda"=lambda, "C"=C, "poly"=poly_order,
                       "auc"=mean_auc, "mad(auc)"=mad_auc, "ams"=mean_ams, "mad(ams)"=mad_ams)
  return(result)
}

# set parameter combinations
logspace <- function(min, max, n) {
  exp(seq(log(min), log(max), length.out = n))
}

n_rbf <- c(0, 2, 4)
lambda <- logspace(1e-4, 20, 3)
C <- c(NA, 1)
poly_order <- 1:3

params <- expand.grid(n_rbf, lambda, C, poly_order)
colnames(params) <- c("n_rbf", "lambda", "C", "poly_order")

# run them all
results <- NULL
for(i in 1:nrow(params)){
  print(i)
  res <- run_models(K=5, n_rbf = params[i,1], lambda = params[i,2], C=params[i,3], poly_order=params[i,4])
  results <- rbind(results, res)
}

write.csv(results, file="output/results_6.csv")
