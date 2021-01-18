
## ---- echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(fs)
library(purrr)

script_dir <- dirname(sys.frame(1)$ofile)
# if (!path_has_parent(script_dir, "R")) {
#   script_dir <- path_join(c(script_dir, "R"))
#   } else {
#     script_dir <- substr(script_dir, 1, nchar(script_dir)-2)
# }
setwd(script_dir)

source("utility_funcs.r")
source("plot_funcs.r")
source("lr_funcs.r")
source("lr1_funcs.r")

## ----load data---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# divide data into training kaggle set, and retain hold-out (before further cross-validation partitioning)
source_path <- path_join(c(path_dir(path_dir(getwd()))))
filename <- "atlas-higgs-challenge-2014-v2.csv"
filepath <- path_join(c(source_path, "LHC_dump", "R", filename))
data <- import_data(filepath)

run_models <- function(model_init, data, train_label, val_label, K, G, n_rbf, lambda, C=1, results_label=2, iter=0, thresholds=c(0.4,0.4,0.4), do_poly_transform=FALSE) {

  ## ----training set------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  train_idx <- get_subset_idx(data$kaggle_s, train_label)
  X <- data$X[train_idx, ]
  y <- data$y[train_idx]
  # need to use kaggle weignts to make AMS correct?!
  kaggle_w <- data$kaggle_w[train_idx]
  w <- data$kaggle_w[train_idx]
  nj <- data$nj[train_idx]
  e_id <- data$e_id[train_idx]


  ## ----validation set----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # public leaderboard set
  val_idx <- get_subset_idx(data$kaggle_s, val_label)
  Xv <- data$X[val_idx, ]
  yv<- data$y[val_idx]
  wv <- data$kaggle_w[val_idx]
  njv <- data$nj[val_idx]


  ## ----add features------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # modify features
  X <- reduce_features(X)
  X <- invert_angle_sign(X)
  Xv <- reduce_features(Xv)
  Xv <- invert_angle_sign(Xv)

if (do_poly_transform) {
  X <- poly_transform(X, 2)
  Xv <- poly_transform(Xv, 2)
}

  # ensure X and Xv have the same columns
  cols2keep <- setdiff(colnames(Xv), setdiff(colnames(Xv), colnames(X)))
  Xv <- Xv[, cols2keep]

  Xv <- scale_dat(Xv, X, na.rm=TRUE)


  ## ----init params-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  s <- avg_median_pairwise_distance(X)

  # K-Fold CV partitioning
  kI <- partition_data(length(y), K, random = TRUE)

  # number of rows and columns of data matrix
  n <- nrow(X)
  d <- ncol(X) + n_rbf

  # sum weights for to renormalise for AMS in partitions
  sum_w <- sum(w)

  # set colours for jet groups
  colours <- generate_colours(G)

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
  b <- matrix(, nrow=d, ncol=G*K)

  ams <- rep(NA, length=G*K)
  auc <- rep(NA, length=G*K)

  par(mfrow=c(2, 3))
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

      col_idx <- get_valid_cols(colnames(X), features_to_rm, mj)

      Xtrain <- as.matrix(Xtrain[, col_idx])
      Xtest <- as.matrix(Xtest[, col_idx])

      Xtest <- scale_dat(Xtest, Xtrain)
      Xtrain <- scale_dat(Xtrain, Xtrain)

      model_idx <- get_model_idx(mj, k, K)

      #fit a logistic regression model to the CV training data
      models[[model_idx]] <- model_init(X=Xtrain, y=y[fit_row_idx])

      #use it to predict the classifications of the test data
      p_hat <- models[[model_idx]]$predict(Xtest)

      #create an ROC curve object
      rocs[[model_idx]] <- ROC_curve$new(y[test_row_idx], p_hat)

      ams_obj[[model_idx]] <- AMS_data$new(y[test_row_idx], p_hat, w[test_row_idx], sum_w=sum_w)
    }
  }


  ## ----avg auc and plot roc, error=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------
  auc <- sapply(rocs, function(x) x$calc_auc())

  ## ----avg ams and plot roc, error=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------
  ams <- sapply(ams_obj, function(x) x$calc_ams())

  ## ----full model--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # loop over folds
  #create lists to hold the k models and k roc curves
  full_models <- vector("list", G)
  full_rocs <- vector("list", G)
  rbf_idx <- matrix(, nrow=G, ncol=n_rbf)
  # check warnings?
  for (mj in 1:G) {
    # loop over sets of jet number {0, 1, 2+} and mH presence/absence
    j <- jet_cats[mj]
    fit_row_idx <- idx_jet_cat(nj, j) & idx_higgs_mass(X, mj, G)

    if (n_rbf > 0) {
      # add r RBF centroid features, using the same reference centroids in training and testing sets
      rbf_centroids <- get_rbf_centroids(X[fit_row_idx, ], n_rbf)
      Xi <- rbf_centroids$"xi"
      # record which rows to use for OOS in public set (or any other OOS)
      rbf_idx[mj, ] <- rbf_centroids$"idx"
      Xtrain <- add_rbf_features(X[fit_row_idx, ], s, n_rbf, Xi=Xi)
    } else {
      Xtrain <- X[fit_row_idx, ]
    }

    col_idx <- get_valid_cols(colnames(Xtrain), features_to_rm, mj)

    Xtrain <- as.matrix(Xtrain[, col_idx])
    Xtrain <- scale_dat(Xtrain, Xtrain)

    model_idx <- get_model_idx(mj, 1, 1)

    #fit a logistic regression model to the CV training data
    full_models[[model_idx]] <- model_init(X=Xtrain, y=y[fit_row_idx])
  }


  ## ----validation--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  amsv_obj <- vector("list", G)

  sum_wv <- sum(wv)
  for (mj in 1:G) {
    j <- jet_cats[mj]
    test_row_idx <- idx_jet_cat(njv, j) & idx_higgs_mass(Xv, mj, G)

    if (n_rbf > 0) {
      # add r RBF centroid features, using the same reference centroids in training and testing sets
      Xi <- X[rbf_idx[mj, ], ]
      Xtest <- add_rbf_features(Xv[test_row_idx, ], s, n_rbf, Xi=Xi)
    } else {
      Xtest <- Xv[test_row_idx, ]
    }

    col_idx <- get_valid_cols(colnames(Xtest), features_to_rm, mj)

    Xtest <- as.matrix(Xtest[, col_idx])

    model_idx <- get_model_idx(mj, 1, 1)

    #use it to predict the classifications of the test data
    p_hat <- full_models[[model_idx]]$predict(Xtest)

    #create an ROC curve object
    full_rocs[[model_idx]] <- ROC_curve$new(yv[test_row_idx], p_hat)

    amsv_obj[[model_idx]] <- AMS_data$new(yv[test_row_idx], p_hat, wv[test_row_idx], sum_w=sum_wv)
  }


  ## ----plot validation auc-----------------------------------------------------------------------------------------------------------------------------------------------------------------
  . <- sapply(full_rocs, function(x) x$calc_auc())
  aucv <- sapply(full_rocs, function(x) round(x$auc, 3))

  ## ----plot validation ams thresholds------------------------------------------------------------------------------------------------------------------------------------------------------
  amsv <- sapply(amsv_obj, function(x) x$calc_ams())

  ## ----print-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  sprintf("Results for lambda=%.2g, G=%i, n_rbf=%i, K=%i on training set ('%s') and validation set ('%s')", lambda, G, n_rbf, K, train_label, val_label)
  sprintf("CV OOS AUC = %.2f ± %.1f", mean(auc), mad(auc))
  sprintf("CV OOS AMS = %.2f ± %.1f", mean(ams), mad(ams))
  sprintf("Validation set AUC = %.2f ± %.1f", mean(aucv), mad(aucv))
  sprintf("Validation set AMS = %.2f ± %.1f", mean(amsv), mad(amsv))

  ## ----record to file----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  filename <- sprintf("results_%s.csv", results_label)
  model_type <- class(models[[1]])[1]
  if (!file.exists(filename)){
    result <- list("Train"=train_label, "Validation"=val_label, "K"=K, "G"=G, "model_type"=model_type, "n_rbf"=n_rbf, "lambda"=lambda, "c"=C, "poly"=do_poly_transform, "auc"=mean(auc), "mad(auc)"=mad(auc), "ams"=mean(ams), "mad(ams)"=mad(ams), "aucv"=mean(aucv), "mad(aucv)"=mad(aucv), "amsv"=mean(amsv), "mad(amsv)"=mad(amsv))
    write.table(as.data.frame(result), file = filename, append = TRUE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = TRUE, qmethod = c("escape", "double"))
  } else {
    result <- list("Train"=train_label, "Validation"=val_label, "K"=K, "G"=G, "model_type"=model_type, "n_rbf"=n_rbf, "lambda"=lambda, "c" = C, "poly"=do_poly_transform, "auc"=mean(auc), "mad(auc)"=mad(auc), "ams"=mean(ams), "mad(ams)"=mad(ams), "aucv"=mean(aucv), "mad(aucv)"=mad(aucv), "amsv"=mean(amsv), "mad(amsv)"=mad(amsv))
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
lambda <- logspace(1e-4, 100, 5)
C <- logspace(1, 100, 5)
thresholds <- c(0.4, 0.4, 0.4)
do_poly_transform = FALSE
# L1 logistic regression (using CVXR)
# rm(logistic)
# model_init <- partial(logistic_l1_model$new, C=C)
results_label <- "experiments1"

library(foreach)
library(doParallel)
registerDoParallel(3)

n_grid_2 <- length(C)

# run them all
foreach (i=1:length(n_rbf)) %do% {
  for (ii in 1:n_grid_2) {
    print(10*(i-1) + ii)
    # model_init <- partial(logistic_model$new, lambda=lambda[ii])
    model_init <- partial(logistic_l1_model$new, C=C[ii])
    run_models(model_init, data, train_label, val_label, K, G, n_rbf[i], lambda[ii], C[ii], results_label=results_label, iter=10*(i-1) + ii, thresholds=thresholds, do_poly_transform=do_poly_transform)
  }
}
