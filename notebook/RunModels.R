library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(fs)
library(purrr)
library(lhc)

## ----load data---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# divide data into training kaggle set, and retain hold-out (before further cross-validation partitioning)
source_path <- path_join(c(path_dir(getwd())))
filename <- "atlas-higgs-challenge-2014-v2.csv"
filepath <- path_join(c(source_path, "LHC_dump", "R", filename))
data <- import_data(filepath)

# function just to loop over fitting procedure with various input parameters to faciliate running experiments over parameter space (e.g. a gridsearch)
run_models <- function(model_init, data, train_label, val_label, K, G, n_rbf, lambda, C=1, results_label=2, poly_order=1) {

  ## ----training set------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  train_idx <- get_subset_idx(data$kaggle_s, train_label)
  X <- data$X[train_idx, ]
  y <- data$y[train_idx]
  # need to use kaggle weignts to make AMS correct?!
  kaggle_w <- data$kaggle_w[train_idx]
  w <- data$kaggle_w[train_idx]
  nj <- data$nj[train_idx]
  e_id <- data$e_id[train_idx]


  ## ----add features------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # modify features
  X <- reduce_features(X)
  X <- invert_angle_sign(X)

  # add polynomial transformations of order poly_order to our covariates
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

  # sum weights for to renormalise for AMS in partitions
  sum_w <- sum(w)

  ## ----jet/missing-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # get missing rows. separate by number of jets and presence of Higgs mass and fit separate models
  # find columns with features with any missing values for each number of jets: 0, 1, 2+ in combination with the presence (or absence) of the Higgs mass, defined by j=1,2,3, 4, 5, 6. i.e. j=1 => nj=0 & mH != -999, j=2 =>
  jet_cats <- c(1:3, 1:3)
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
      fit_row_idx <- kI != k & idx_jet_cat(nj, j) & idx_higgs_mass(X, mj, G)
      test_row_idx <- kI == k & idx_jet_cat(nj, j) & idx_higgs_mass(X, mj, G)

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
      models[[model_idx]] <- model_init(X=Xtrain, y=y[fit_row_idx])

      #use it to predict the classifications of the test data
      p_hat <- models[[model_idx]]$predict(Xtest)

      #create an ROC curve object
      rocs[[model_idx]] <- ROC_curve$new(y[test_row_idx], p_hat)

      # create an AMS object
      w_this_partition <- w[test_row_idx] * sum_w/sum(w[test_row_idx])
      ams_obj[[model_idx]] <- AMS_data$new(y[test_row_idx], p_hat, w_this_partition)
    }
  }


  ## ----avg auc and plot roc----------------------------------------------------------------------------------------------------------------------------------------------------
  auc <- sapply(rocs, function(x) x$calc_auc())

  ## ----avg ams and plot roc----------------------------------------------------------------------------------------------------------------------------------------------------
  ams <- sapply(ams_obj, function(x) x$calc_ams())

  get_mean_mad_ams <- function(ams, G, K) {
    mads <- rep(NA, G)
    for (g in 1:G) {
      model_set <- ((g-1)*K+1):(g*K)
      mads[g] <- mean(apply(ams[, model_set], 1, mad))
    }
    return(mean(mads))
  }

  mad_ams <- get_mean_mad_ams(ams, G, K)

  ## ----record to file----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  filename <- sprintf(path_join(c(getwd(), "output/results_%s.csv")), results_label)
  model_type <- class(models[[1]])[1]
  if (!file.exists(filename)){
    result <- list("Train"=train_label, "Validation"=val_label, "K"=K, "G"=G, "model_type"=model_type, "n_rbf"=n_rbf, "lambda"=lambda, "c"=C, "poly"=poly_order, "auc"=mean(auc), "mad(auc)"=mad(auc), "ams"=mean(ams), "mad(ams)"=mad_ams)#, "aucv"=mean(aucv), "mad(aucv)"=mad(aucv), "amsv"=mean(amsv), "mad(amsv)"=mad(amsv))
    write.table(as.data.frame(result), file = filename, append = TRUE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = TRUE, qmethod = c("escape", "double"))
  } else {
    result <- list("Train"=train_label, "Validation"=val_label, "K"=K, "G"=G, "model_type"=model_type, "n_rbf"=n_rbf, "lambda"=lambda, "c" = C, "poly"=poly_order, "auc"=mean(auc), "mad(auc)"=mad(auc), "ams"=mean(ams), "mad(ams)"=mad_ams)#, "aucv"=mean(aucv), "mad(aucv)"=mad(aucv), "amsv"=mean(amsv), "mad(amsv)"=mad(amsv))
    write.table(as.data.frame(result), file = filename, append = TRUE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"))
  }
}

# set parameter combinations
logspace <- function(min, max, n) {
  exp(seq(log(min), log(max), length.out = n))
}

train_label <- c("t")
val_label <- c("b")
K <- 10
G <- 3
n_rbf <- 0:5
lambda <- logspace(1e-4, 20, 15)
C <- 1#logspace(1, 100, 5)
poly_order <- 3
# L1 logistic regression (using CVXR)
# rm(logistic)
# model_init <- partial(logistic_l1_model$new, C=C)
results_label <- "experiments5"

library(foreach)
library(doParallel)
registerDoParallel(4)

n_grid_2 <- length(lambda)

# run them all
for (poly_order in 1:3) {
  foreach (i=1:length(n_rbf)) %dopar% {
    for (ii in 1:n_grid_2) {
      print(c(poly_order, i, ii))
      model_init <- partial(logistic_model$new, lambda=lambda[ii])
      # model_init <- partial(logistic_l1_model$new, C=C[ii])
      run_models(model_init, data, train_label, val_label, K, G, n_rbf[i], lambda[ii], C[ii], results_label=results_label, poly_order=poly_order)
    }
  }
}
