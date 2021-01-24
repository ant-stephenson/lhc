#' Import raw LHC data
#'
#' This function provides a standard way to load in the LHC dataset for our analysis pipeline.
#' Undefined data (-999s) are replaces with NAs, columns are split into variables,
#' labels and supplementary data.
#'
#' @param filepath str location of csv
#' @return named list of X, y, w, kaggle_w, kaggle_s, e_id, nj
#' @importFrom data.table fread
#' @export
import_data <- function(filepath="atlas-higgs-challenge-2014-v2.csv") {
  # Fast load large csv with fread from data.table
  lhc_data <- as.data.frame(fread(filepath))

  # Assign target variable
  y <- rep(0, nrow(lhc_data))
  y[lhc_data$Label == "s"] <- 1

  # Assign other useful label vectors (non-variable data)
  w <- lhc_data$Weight
  kaggle_w <- lhc_data$KaggleWeight
  kaggle_s <- lhc_data$KaggleSet
  e_id <- lhc_data$EventId
  nj <- lhc_data$PRI_jet_num

  # Setup design matrix
  X_header <- names(lhc_data)[grep("Weight|Label|KaggleWeight|KaggleSet|EventId|PRI_jet_num", names(lhc_data), invert=TRUE)]
  X <- lhc_data[, X_header]

  # Replace -999s with NAs
  X[X==-999] <- NA

  output <- list(X=X, y=y, w=w, kaggle_w=kaggle_w, kaggle_s=kaggle_s, e_id=e_id, nj=nj)
  return(output)
}

#' Get boolean vector of rows with j=0,1 or 2+
#' @param nj Vector of number of jets for each point
#' @param j jet group
#' @export
idx_jet_cat <- function(nj, j) {
  if (j == 1) {
    nj == 0
  } else if (j == 2) {
    nj == 1
  } else if (j == 3) {
    nj >= 2
  } else {stop("Incorrect jet category specified")}
}

#' Get boolean index for rows with missing/or not missing (depending on G/j) Higgs mass
#' @param X matrix of covariates
#' @param j jet group
#' @param G number of jet groups
#' @return vector of bools
#' @export
idx_higgs_mass <- function(X, j, G) {
  # for j>3 take rows with Higgs missing (and drop Higgs)
  if (G == 6) {
    is_missing <- is.na(X$"DER_mass_MMC")
    if (j > 3) {
      return(!is_missing)
    } else {return(is_missing)}
  } else {
    return(rep(TRUE, dim(X)[1]))
  }
}

#' create list of feature names we want to omit based on jet group, or constant values/missing
#' @param X matrix of covariates
#' @param G number of jet groups
#' @param kI fold indices (test label)
#' @param nj Vector of number of jets for each point
#' @return nested list of column names
#' @export
set_features_to_rm <- function(X, G, kI, nj) {
  jet_cats <- rep(c(1:3), G/3)
  K <- max(kI)
  features_to_rm <- list(list(), list(), list())
  for (mj in 1:G) {
    j <- jet_cats[mj]
    .features_to_rm <- colnames(X)[colSums(is.na(X[idx_jet_cat(nj, j) & idx_higgs_mass(X, mj, G),])) > 0]

    # check if any features are constant over a fold of a jet group as then we'll have multicollinearity issues
    const_features <- list()
    for (k in 1:K) {
      fit_row_idx <- kI != k & idx_jet_cat(nj, j) & idx_higgs_mass(X, mj, G)
      const_features <- union(const_features, get_const_features(X[fit_row_idx,]))
    }

    features_to_rm[[mj]] <- union(.features_to_rm, unlist(const_features))
  }
  return(features_to_rm)
}


#' Reduce feature space dimensionality by exploiting redundancy
#'
#' @param X matrix of covariates
#' @return X augmented matrix of covariates
#' @export
# based on
# https://www.kaggle.com/c/higgs-boson/discussion/9576
reduce_features <- function(X) {
  X$"PRI_lep_phi-PRI_tau_phi" <- (X[, "PRI_lep_phi"] - X[, "PRI_tau_phi"]) %% (2*pi)
  X$"PRI_met_phi-PRI_tau_phi" <- (X[, "PRI_met_phi"] - X[, "PRI_tau_phi"]) %% (2*pi)
  X$"PRI_jet_leading_phi-PRI_tau_phi" <- (X[, "PRI_jet_leading_phi"] - X[, "PRI_tau_phi"]) %% (2*pi)
  X$"PRI_jet_subleading_phi-PRI_tau_phi" <- (X[, "PRI_jet_subleading_phi"] - X[, "PRI_tau_phi"]) %% (2*pi)

  to_drop <- c("PRI_tau_phi", "PRI_lep_phi", "PRI_met_phi", "PRI_jet_leading_phi", "PRI_jet_subleading_phi")
  return(X[ , !(names(X) %in% to_drop)])

  # new_phi=(rot_phi+3*pi) %% (2*pi) - pi
}

#' Invert Angle Sign
#'
#' Uses the sign of the pseudorapidity of the tau particle to modify the sign of the pseudorapidity of the leptons and jets
#' on the basis that the interaction should be invariant to rotations of pi about the beam (z) axis.
#' \deqn{\eta(\theta) = -\log \tan \frac{\theta}{2}}
#' \deqn{\eta(\pi-\theta) = -\eta(\theta)}
#' @export
invert_angle_sign <- function(X) {
  signs <- sign(X$"PRI_tau_eta")
  X$"PRI_tau_eta" <- signs * X$"PRI_tau_eta"
  X$"PRI_lep_eta" <- signs * X$"PRI_lep_eta"
  X$"PRI_jet_leading_eta" <- signs * X$"PRI_jet_leading_eta"
  X$"PRI_jet_subleading_eta" <- signs * X$"PRI_jet_subleading_eta"
  return(X)
}
