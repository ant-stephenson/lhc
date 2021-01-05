#' Select indices of rows in x which correspond to values of labels, which can be multiple elements
#' @param x vector of labels
#' @param labels reference labels to compare to
#' e.g. x <- c("t", "b", "t", "v", "b"); labels <- c("t", "b")
#' get_subset_idx(x, labels) = T, T, T, F, T
get_subset_idx <- function(x, labels) {
  !is.na(match(x, labels))
}

#' Import data from the Higgs csv file. Reorganise into covariates, labels and supplementary data (weights, ids)
#' Replace -999s with NA
#' Standardize covariates
#' 
#' @param filepath str location of csv
#' @return named list of X, y, w, kaggle_w, kaggle_s, e_id, nj
import_data <- function(filepath="atlas-higgs-challenge-2014-v2.csv") {
    # Import data - put all this in a function to create a pipeline
    lhc_data <- read.csv(filepath)

    n <- nrow(lhc_data)

    # Assign target variable
    label <- lhc_data$Label
    y <- rep(0,n)
    y[label == "s"] <- 1

    # Assign other useful label vectors (non-variable data)
    w <- lhc_data$Weight
    kaggle_w <- lhc_data$KaggleWeight
    kaggle_s <- lhc_data$KaggleSet
    e_id <- lhc_data$EventId
    nj <- lhc_data$PRI_jet_num

    # Setup design matrix, with intercept, taking care with -999s: this is why they're annoying
    X_header <- names(lhc_data)[grep("Weight|Label|KaggleWeight|KaggleSet|EventId|PRI_jet_num", names(lhc_data), invert=TRUE)]
    X <- as.matrix(lhc_data[, X_header])

    # shouldn't be necessary, but record -999s locations. makes it easy to transform back to -999 if desired. for now, replace with NA which should propagate
    loc999s <- X == -999
    X[loc999s] <- NA

    X <- data.frame(X)

    output <- list(X=X, y=y, w=w, kaggle_w=kaggle_w, kaggle_s=kaggle_s, e_id=e_id, nj=nj)

    return(output)
}

#define a function to scale features of a matrix with reference to another matrix
#useful because you can normalise X_train, and apply the same transformation to X_test
#not designed for data with -999s! 
scale_dat <- function(X, ref, na.rm=FALSE){
  if(ncol(X) != ncol(ref)) stop('Two inputs must have the same number of columns')

  #calculate column means and sds of ref, ignoring NAs
  mu <- colMeans(ref, na.rm=na.rm)
  sd <- apply(ref, 2, function(x) sd(x, na.rm=na.rm))
  
  #transform columns of X
  for(i in 1:ncol(ref)){
    X[,i] <- (X[,i] - mu[i]) / sd[i] #is there a smarter way to do this not in a loop?
  }

  #also add column of 1s called intercept
  Intercept <- rep(1, nrow(X))
  X <- cbind(Intercept, X)
  return(X)
}

#' Get boolean vector of rows with j=0,1 or 2+
#' @param nj Vector of number of jets for each point
#' @parma j jet group
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
#' @parma j jet group
#' @param G number of jet groups
#' @return vector of bools
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
set_features_to_rm <- function(X, G, kI, nj) {
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

#' Uses the sign of the pseudorapidity of the tau particle to modify the sign of the pseudorapidity of the leptons and jets, on the basis that the interaction should be invariant to rotations of $\pi$ about the beam (z) axis.
#' $\eta(\theta) = -\log\tan\frac{\theta}{2}$
#' $\eta(\pi-\theta) = -\eta(\theta)$
invert_angle_sign <- function(X) {
    signs <- sign(X$"PRI_tau_eta")
    X$"PRI_tau_eta" <- signs * X$"PRI_tau_eta"
    X$"PRI_lep_eta" <- signs * X$"PRI_lep_eta"
    X$"PRI_jet_leading_eta" <- signs * X$"PRI_jet_leading_eta"
    X$"PRI_jet_subleading_eta" <- signs * X$"PRI_jet_subleading_eta"
    return(X)
}

#' Add n_centroids RBF features to X
#' 
#' @param X covariate matrix
#' @param s median pairwise distance of points in X
#' @param n_centroids number of RBF features to add
#' @param Xi matrix of reference points to calculate RBF w.r.t 
#' @return X augmented covariate matrix
add_rbf_features <- function(X, s, n_centroids, Xi=NULL) {
    set.seed(1234)
    if (is.null(Xi)) {
        rbf_centroids <- get_rbf_centroids(X, n_centroids)
        Xi <- rbf_centroids$"xi"
    }
    for (i in 1:n_centroids) {
        label <- sprintf("rbf%i", i)
        X[, label] <- rbf_feature(X, s, xi=Xi[i,])
    }
    return(X)
}

#' Get reference points for RBF centroids
#' 
#' @param X covariate matrix
#' @param n_centroids number of RBF centroids
#' @param idx [Optional] location of RBF centroid reference points
#' @return list of Xi (centroid points) and idx (location of them)
get_rbf_centroids <- function(X, n_centroids, idx=NULL) {
    n <- nrow(X)
    d <- ncol(X)
    if (is.null(idx)) {
        idx <- sample(1:n, n_centroids)
    }
    Xi <- matrix(, nrow=n_centroids, ncol=d)
    for (i in 1:n_centroids) {
        Xi[i, ] <- unlist(X[idx[i],])
    }
    return(list("xi"=Xi, "idx"=idx))
}

#' Compute single RBF feature at some centroid i in idx (or xi in Xi) 
#' 
#' @param X covariate matrix
#' @param s median pairwise distance of points in X
#' @param idx [Optional] location of reference centroid
#' @param xi [Optional] reference centroid
rbf_feature <- function(X, s, idx=NULL, xi=NULL) {
    n <- nrow(X)
    if (is.null(idx)) {
        set.seed(1234)
        idx <- sample(1:n, 1)
    }
    if (is.null(xi)) {
        xi <- X[idx, ]
    }
    Xi <- matrix(xi, nrow=n, ncol=length(xi), byrow=TRUE)
    Xi <- apply(Xi, 2, unlist)
    rbf <- exp(-pairwise_distance(X, Xi)/(2 * s))
    return(rbf)
}

#' Calculate average of median pairwise distances for between all adjacent points
#' 
#' @param X covariate matrix
#' @return s median pairwise distance
avg_median_pairwise_distance <- function(X) {
    n <- nrow(X)
    X_perm <- X
    s <- rep(0,2)
    for (i in 1:2) {
        X_perm <- permute_matrix(X_perm)
        s[i] <- median(pairwise_distance(X, X_perm))
    }
    mean(s)
}

#' Calculate distance for each row of X0 and X1
#' 
#' @param X0 covariate matrix
#' @param X1 covariate matrix
#' @return vector of distances
pairwise_distance <- function(X0, X1) {
    rowSums((X0 - X1)^2, na.rm=TRUE)
}

#' Cyclic permutation of rows of X by r rows
#' 
#' @param X covariate matrix
#' @param r number of rows to permute (default=1)
permute_matrix <- function(X, r=1) {
    n <- nrow(X)
    X_perm <- rbind(X[(n-r+1):n,], X[1:(n-r),])
    return(X_perm)
}

poly_transform <- function(X, b=2){
    for(i in 2:b){
        Xb <- apply(X[, ], 2, function(col) col^b)
        colnames(Xb) <- paste0(colnames(Xb), "^", b)
        X <- cbind(X, Xb)
    }
    #remove highly correlated variables
    cors <- cor(X)
    cors[!lower.tri(cors)] <- 0
    X <- X[, !apply(cors,2,function(x) any(abs(x) > 0.80, na.rm=TRUE))]
    return(X)
}

#' find colnames for columns that are constant (e.g. all 1, -999, NA etc)
#' @param X matrix of covariates
#' @return list of column names
get_const_features <- function(X) {
    col_sd <- apply(X, 2, sd)
    cols <- col_sd == 0 | is.na(col_sd)
    return(setdiff(colnames(X)[cols], "Intercept"))
}

#' Partition data into (random) folds for cross-validation. 
#' 
#' @param n number of rows
#' @param k number of folds
#' @param random flag to choose whether to randomly select
#' @return ind vector of integers denoting the OOS fold each row belongs to
#' Returns vector of indices denoting the OOS index, i_e. for rows with I_i=1, those are OOS for i=1
partition_data <- function(n, k, random=FALSE){
    #minimum size of the group
    size <- floor(n/k)
    
    #number of remaining points after even division
    remainder <- n %% k
    
    #repeat 1:k size times, then 1:remainder (makes a vector length n)
    ind <- c(rep(seq(1:k), size), seq(1, remainder, length.out = remainder))
    
    #if you want, shuffle the index
    if(random==T){
        ind <- ind[sample(ind)]
    }
}

#' helper function to get the index corresponding to the model built on folds {l !=k} and for jet number j {1,2,3}, ordering columns by fold and then with nesting on j i.e. first six cols are k=1, j=1,2,3; k=2, j=1,2,3; etc
#' @param j this jet group
#' @param k this fold
#' @param K number of folds
get_model_idx <- function(j, k, K) {
    K*(j-1) + k
}

#' does the inverse procedure to get_model_idx, s.t. if idx <- get_model_idx(j, k, K)
#' then (j,k) <- inv_model_idx(idx) (if R output tuples)
#' @param idx model index
#' @param K number of folds
#' @return numeric pair of j and k
inv_model_idx <- function(idx, K) {
    # assumes indexing loops internally over j and externally over k
    j <-  ceiling(idx/K)
    k <- idx - K*(j-1)
    return(c(j, k))
}

#' Gets the columns from header that aren't in features_to_rm
#' @param header list of column names e.g. names(X)
#' @param features_to_rm list of feature names we want to remove
#' @param j jet group
#' @return index of columns to retain
get_valid_cols <- function(header, features_to_rm, j) {
    match(setdiff(header, features_to_rm[[j]]), header)
}

#' Calculate logistic function
#' @param x float
#' @return logistic(x)
logistic <- function(x) {
    1/(1 + exp(-x))
}

#' Calculate logit function
#' @param p float in [0,1]
#' @return logit(p)
logit <- function(p) {
    log(p/(1-p))
}

#' Define polynomial kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @param b order of polynomial
#' @return scalar of (1+x^Tx)^b
poly_kernel <- function(x_i, x_j, b) {
    return((t(x_i) %*% x_j + 1)^b)
}

#' Define trigonometric kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @param b order of polynomial
trig_kernel <- function(x_i, x_j, b=0) {
    output <- 0
    for (k in 0:b) {
        output <- output + cos(k * (x_i - x_j))
    }
    return(sum(output))
}

#' Define linear kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
lin_kernel <- function(x_i, x_j, ...) {
    return(t(x_i) %*% x_j)
}

#' RBF kernel
#' @param x_i point in R^d
#' @param x_j point in R^d
#' @param sigma bandwidth hyperparameter
rbf_kernel <- function(x_i, x_j, sigma) {
    return(exp(-t(x_i - x_j) %*% (x_i - x_j) / (2 * sigma^2)))
}

#' Function factory to partially call kernel function to return the kernel function with it's hyperparameters set
#' @param ckernel kernel function with args (x_i,x_j,hyper)
#' @return function with args (x_i,x_j)
tuned_kernel <- function(ckernel, ...) {
    function(x_i, x_j){
        ckernel(x_i, x_j, ...)
    }
}

#' Compute the kernel matrix over the (training) set X
#' @param X covariate matrix (nxd)
#' @param ckernel kernel function with args (x_i, x_j)
#' @return K kernel matrix (nxn)
calc_K <- function(X, ckernel){
    n <- nrow(X)
    K <- matrix(, nrow=n, ncol=n)
    for (i in 1:n){
        for (j in 1:n){
            K[i,j] <- ckernel(X[i,], X[j,])
        }
    }
    return(K)
}

#' Calculate predictions for test points by computing the kernel matrix over the training points and the kernel vector(/matrix) over test and training points
#' 
#' @param X_test matrix of covariate test points
#' @param X_train matrix of covariate training points
#' @param y response vector
#' @param L ?
#' @param ckernel kernel function k(x_i,x_j)
svm_predict <- function(X_test, X_train, y, L, ckernel=lin_kernel) {
    X_test <- as.matrix(X_test)
    X_train <- as.matrix(X_train)
    n_test <- nrow(X_test)
    n_train <- nrow(X_train)
    k <- matrix(, nrow=n_train, ncol=n_test)
    K <- matrix(, nrow=n_train, ncol=n_train)
    for (i in 1:n_train){
        for (m in 1:n_test){
            k[i,m] <- ckernel(X_test[m,], X_train[i,])
        }
        for (j in 1:n_train){
            K[i,j] <- ckernel(X_train[i,], X_train[j,])
        }
    }
    f <- t(k) %*% solve(K + L * diag(n_train)) %*% y
    return(f)
}

#' Thresholding function
#' 
#' @param p vector of probabilities
#' @param thresh threshold over which we assign output of 1
decide <- function(p, thresh = 0.5) {
    label <- as.numeric(p > thresh)
    return(label)
}

#' the AMS metric. 
#' note s = sum_{i in B \cup G}w_i and b = sum_{i in B \cup G}w_i; 
#' i.e. the sum of the weights of succesful signal classifications (TP) 
#' and the sum of the weights of incorrect signal classifications (FP) respectively
#' @param s count of true positives
#' @param b count of false positives
ams_metric <- function(s, b) {
    br <- 10
    sqrt(2 * ((s + b + br) * log(1 + s/(b + br)) - s))
}

#' Count the number of true positives
#' need to make sure dims of y, y_hat and w are the same
#' @param y response vector
#' @param y_hat predicted response vector
#' @param w weights
count_s <- function(y, y_hat, w) {
    stopifnot(length(y) == length(y_hat))
    stopifnot(length(y) == length(w))
    s <- sum(w[y == 1 & y_hat == 1])
    return(s)
}

#' Count the number of false positives
#' need to make sure dims of y, y_hat and w are the same
#' @param y response vector
#' @param y_hat predicted response vector
#' @param w weights
count_b <- function(y, y_hat, w) {
    stopifnot(length(y) == length(y_hat))
    stopifnot(length(y) == length(w))
    b <- sum(w[y == 0 & y_hat == 1])
    return(b)
}

#' Calculate the AMS
#' @param y response vector
#' @param y_hat predicted response vector
#' @param w weights
#' @param sum_w total sum of weights for renormalisation
#' @return ams
calculate_ams_partition <- function(y, y_hat, w, sum_w=1) {
    w <- w * sum_w/sum(w)
    y_hat <- as.numeric(y_hat)
    s <- count_s(y, y_hat, w)
    b <- count_b(y, y_hat, w)
    ams <- ams_metric(s, b)
    return(ams)
}


