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

    # standardize data and add intercept - not technically OOS if we do this on all data
    X <- scale(X, TRUE, TRUE)
    n <- nrow(X)
    X <- cbind(rep(1, n), X)
    X_header <- c("Intercept", X_header)
    X <- data.frame(X)

    output <- list(X=X, y=y, w=w, kaggle_w=kaggle_w, kaggle_s=kaggle_s, e_id=e_id, nj=nj)

    return(output)
}

# reduce feature space
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

invert_angle_sign <- function(X) {
    signs <- sign(X$"PRI_tau_eta")
    X$"PRI_tau_eta" <- signs * X$"PRI_tau_eta"
    X$"PRI_lep_eta" <- signs * X$"PRI_lep_eta"
    X$"PRI_jet_leading_eta" <- signs * X$"PRI_jet_leading_eta"
    X$"PRI_jet_subleading_eta" <- signs * X$"PRI_jet_subleading_eta"
    return(X)
}

# add cake variables. doesn't work at the moment due to no overlapping event ids?!
add_cake_features <- function(X, e_id) {
    cake <- read.csv(file = "cake_features.csv")
    X$"EventId" <- e_id
    X <- merge(X, cake[,1:3])
    return(X[, !(names(X) %in% "EventId")])
}

# add n_centroids RBF features to X
add_rbf_features <- function(X, s, n_centroids, Xi=NULL) {
    set.seed(1234)
    if (is.null(Xi)) {
        rbf_centroids <- get_rbf_centroids(X, n_centroids)
        Xi <- rbf_centroids$"xi"
    }
    for (i in 1:n_centroids) {
        label <- sprintf("rbf%i", i)
        X$label <- rbf_feature(X, s, xi=Xi[i,])
    }
    return(X)
}

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

# rbf feature at some centroid i in 
rbf_feature <- function(X, s, idx=NULL, xi=NULL) {
    n <- nrow(X)
    if (is.null(idx)) {
        set.seed(1234)
        idx <- sample(1:n, 1)
    }
    if (is.null(xi)) {
        xi <- X[idx, ]
    }
    Xi <- matrix(rep(as.numeric(xi), n), nrow=n, ncol=ncol(X))
    rbf <- exp(-pairwise_distance(X, Xi)/(2 * s))
    return(rbf)
}

# calculate average of median pairwise distances, cycling over all pairs. this scales _really, really_ badly with n. don't do it!!!
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

# calculate distance for each row of X0 and X1
pairwise_distance <- function(X0, X1) {
    rowSums(apply((X0 - X1)^2, 2, unlist), na.rm=TRUE)
}

# cyclic permutation of rows of X by r rows
permute_matrix <- function(X, r=1) {
    n <- nrow(X)
    X_perm <- rbind(X[(n-r+1):n,], X[1:(n-r),])
    return(X_perm)
}

# Partition data into (random) folds for cross-validation. Returns vector of indices denoting the OOS index, i_e. for rows with I_i=1, those are OOS for i=1
partition_data <- function(X, y, K, random=TRUE) {
    n <- max(nrow(X), ncol(X))
    rows_per_fold <- floor(n/K)
    I <- rep(1, n)
    rem_ind <- 1:n
    for (k in 1:K){
        if (random) {
            ind <- sample(rem_ind, rows_per_fold)
            rem_ind <- setdiff(rem_ind, ind)
        }
        else {
            ind <- ((k - 1) * rows_per_fold + 1):(k * rows_per_fold)
        }
        I[ind] <- k
    }
    return(I)
}

logistic <- function(x) {
    1/(1 + exp(-x))
}

logit <- function(p) {
    log(p/(1-p))
}

# fit a logistic regression model by IRWLS - same thing glm does. probably be useful to make this model into a class with a predict method and then we can interchange it with other models wth that interface
library(Matrix)
logistic_reg <- function(X, y, r=NULL, lambda = 1e-3) {
    invlink <- logistic
    dinvlink <- function(x)  exp(x)/(1 + exp(x))^2
    loss <- function(y, eta) -(y * eta) + log(1 + exp(eta))

    d <- ncol(X)
    n <- nrow(X)

    if (is.null(r)) {r <- rep(1,n)}

    b <- matrix(0, d, 1)
    eta <- matrix(0, n, 1)
    mu <- invlink(eta)
    L <- sum(loss(y, eta)) - lambda * t(b) %*% b

    maxIter <- 10
    tol <- 1e-6

    for (i in 1:maxIter) {
        w <- pmax(dinvlink(eta), 1e-6)
        W <- Diagonal(x = as.numeric(w) * r)
        H <- -t(X) %*% W %*% X + 2 * lambda * diag(rep(1, d))
        grad <- t(X) %*% (y - mu) + 2 * lambda * b

        b <- b - solve(H, grad)
        eta <- X %*% b
        mu <- invlink(eta)
        L_new <- sum(loss(y, eta)) - lambda * t(b) %*% b
        
        if (as.logical(abs(L - L_new) < tol)) {
            print(sprintf("No. iterations: %d", i))
            break
        }
        L <- L_new
    }
    return(b)
}

# think we'll have a problem with data size... nxn will be massive. and it doesn't work
kernel_logistic_reg <- function(K, y, lambda) {
    invlink <- logistic
    dinvlink <- function(x)  exp(x)/(1 + exp(x))^2
    loss <- function(y, eta) -(y * eta) + log(1 + exp(eta))

    n <- nrow(K)

    a <- matrix(0, n, 1)
    eta <- matrix(0, n, 1)
    mu <- invlink(eta)
    L <- sum(loss(y, eta)) - lambda * t(a) %*% K %*% K %*% a

    maxIter <- 10
    tol <- 1e-6

    for (i in 1:maxIter) {
        w <- pmax(dinvlink(eta), 1e-6)
        W <- Diagonal(x = as.numeric(w))
        H <- -K %*% W %*% K + 2 * lambda * K %*% K
        grad <- K %*% (y - mu) + 2 * lambda * K %*% K %*% a

        a <- a - solve(H, grad)
        eta <- K %*% a
        mu <- invlink(eta)
        L_new <- sum(loss(y, eta)) - lambda * t(a) %*% K %*% K %*% a
        
        if (as.logical(abs(L - L_new) < tol)) {
            print(sprintf("No. iterations: %d", i))
            break
        }
        L <- L_new
    }
    return(a)

}

#define ckernels
poly_kernel <- function(x_i, x_j, b) {
    return((t(x_i) %*% x_j + 1)^b)
}

trig_kernel <- function(x_i, x_j, b=0) {
    output <- 0
    for (k in 0:b) {
        output <- output + cos(k * (x_i - x_j))
    }
    return(sum(output))
}

lin_kernel <- function(x_i, x_j, ...) {
    return(t(x_i) %*% x_j)
}

rbf_kernel <- function(x_i, x_j, sigma) {
    return(exp(-t(x_i - x_j) %*% (x_i - x_j) / (2 * sigma^2)))
}

# function to partially call kernel function to return the kernel function with it's hyperparameters set
tuned_kernel <- function(ckernel, ...) {
    function(x_i, x_j){
        ckernel(x_i, x_j, ...)
    }
}

# Compute the kernel matrix over the (training) set X
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

# calculate predictions for test points by computing the kernel matrix over the training points and the kernel vector(/matrix) over test and training points
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

# thresholding function
decide <- function(p, thresh = 0.5) {
    label <- as.numeric(p > thresh)
    return(label)
}

# the AMS metric. note s = sum_{i in B \cup G}w_i and b = sum_{i in B \cup G}w_i; i.e. the sum of the weights of succesful signal classifications (TP) and the sum of the weights of incorrect signal classifications (FP) respectively
ams_metric <- function(s, b) {
    br <- 10
    sqrt(2 * (s + b + br) * log(1 + s/(b + br)) - s)
}

# need to make sure dims of y, y_hat and w are the same
count_s <- function(y, y_hat, w) {
    stopifnot(length(y) == length(y_hat))
    stopifnot(length(y) == length(w))
    s <- sum(w[y == 1 & y_hat == 1])
    return(s)
}

count_b <- function(y, y_hat, w) {
    stopifnot(length(y) == length(y_hat))
    stopifnot(length(y) == length(w))
    b <- sum(w[y == 0 & y_hat == 1])
    return(b)
}

calculate_ams_partition <- function(y, y_hat, w) {
    y_hat <- as.numeric(y_hat)
    s <- count_s(y, y_hat, w)
    b <- count_b(y, y_hat, w)
    ams <- ams_metric(s, b)
    return(ams)
}

calculate_oos_ams <- function(y, y_hat, w, kI) {
    ams <- 0
    for (k in 1:max(kI)) {
        ams <- ams + calculate_ams_partition(y[kI == k], y_hat[kI == k], w[kI == k])
    }
    return(ams)
}

