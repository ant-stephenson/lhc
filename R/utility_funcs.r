#' Select indices of rows in x which correspond to values of labels, which can be multiple elements
#' @param x vector of labels
#' @param labels reference labels to compare to
#' e.g. x <- c("t", "b", "t", "v", "b"); labels <- c("t", "b")
#' get_subset_idx(x, labels) = T, T, T, F, T
#' @export
get_subset_idx <- function(x, labels) {
  !is.na(match(x, labels))
}


#' Add n_centroids RBF features to X
#'
#' @param X covariate matrix
#' @param s median pairwise distance of points in X
#' @param n_centroids number of RBF features to add
#' @param Xi matrix of reference points to calculate RBF w.r.t
#' @return X augmented covariate matrix
#' @export
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
#' @export
get_rbf_centroids <- function(X, n_centroids, idx=NULL) {
    n <- nrow(X)
    d <- ncol(X)
    if (is.null(idx)) {
        idx <- sample(1:n, n_centroids)
    }
    Xi <- matrix(NA, nrow=n_centroids, ncol=d)
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
#' @export
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
#' @export
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
#' @export
pairwise_distance <- function(X0, X1) {
    rowSums((X0 - X1)^2, na.rm=TRUE)
}

#' Cyclic permutation of rows of X by r rows
#'
#' @param X covariate matrix
#' @param r number of rows to permute (default=1)
#' @export
permute_matrix <- function(X, r=1) {
    n <- nrow(X)
    X_perm <- rbind(X[(n-r+1):n,], X[1:(n-r),])
    return(X_perm)
}

#' Run polynomial transform on columns of X (of order b), removing output columns that are highly correlated
#' @param X matrix of covariates
#' @param b order of polynomial
#' @return X augmented matrix of covariates
#' @export
poly_transform <- function(X, b=2){
  orig_cols <- colnames(X)
    for(i in 2:b){
        Xb <- apply(X[, orig_cols], 2, function(col) col^i)
        colnames(Xb) <- paste0(colnames(Xb), "^", i)
        X <- cbind(X, Xb)
    }
    # #remove highly correlated variables - doesnt work
    # cors <- cor(X)
    # cors[!lower.tri(cors)] <- 0
    # X <- X[, !apply(cors,2,function(x) any(abs(x) > 0.80, na.rm=TRUE))]

    #truncate particularly large values that might overflow
    if (b >2) {
      max_val <- apply(X, 2, function(x) log10(max(x, na.rm=TRUE)))
      X <- X[, names(max_val[max_val < 8])]
        # X[abs(X) > 1e8] <- 1e8
    }
    return(X)
}

#' find colnames for columns that are constant (e.g. all 1, -999, NA etc)
#' @param X matrix of covariates
#' @return list of column names
#' @export
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
#' @export
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
#' @export
get_model_idx <- function(j, k, K) {
    K*(j-1) + k
}

#' does the inverse procedure to get_model_idx, s.t. if idx <- get_model_idx(j, k, K)
#' then (j,k) <- inv_model_idx(idx) (if R output tuples)
#' @param idx model index
#' @param K number of folds
#' @return numeric pair of j and k
#' @export
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
#' @export
get_valid_cols <- function(header, features_to_rm, j) {
    match(setdiff(header, features_to_rm[[j]]), header)
}

#' Calculate logistic function
#' @param x float
#' @return logisticf(x)
#' @export
logisticf <- function(x) {
    1/(1 + exp(-x))
}

#' Calculate logit function
#' @param p float in [0,1]
#' @return logit(p)
#' @export
logit <- function(p) {
    log(p/(1-p))
}

#' Thresholding function
#'
#' @param p vector of probabilities
#' @param thresh threshold over which we assign output of 1
#' @export
decide <- function(p, thresh = 0.5) {
    label <- as.numeric(p > thresh)
    return(label)
}

#' Calculate AMS metric
#'
#' \eqn{s = sum_{i in B \cup G}w_i}
#' \eqn{b = sum_{i in B \cup G}w_i};
#' i.e. s is the sum of the weights of successful signal classifications (TP)
#' and b is the sum of the weights of incorrect signal classifications (FP)
#' @param s count of true positives
#' @param b count of false positives
#' @export
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
calculate_ams_partition <- function(y, y_hat, w, sum_w=NULL) {
  if (!is.null(sum_w)) {
    w <- w * sum_w/sum(w)
  }
    y_hat <- as.numeric(y_hat)
    s <- count_s(y, y_hat, w)
    b <- count_b(y, y_hat, w)
    ams <- ams_metric(s, b)
    return(ams)
}


