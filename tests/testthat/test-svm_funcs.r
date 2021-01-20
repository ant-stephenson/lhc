context("Test SVM functions on generated data")

# Create two distinct point clouds
n <- 50
y <- c(rep(1, n), rep(-1, n))
x_1 <- c(rnorm(n, 5, 1), rnorm(n, 10, 1))
x_2 <- c(rnorm(n, 5, 1), rnorm(n, 10, 1))
X <- matrix(c(x_1, x_2), ncol=2)

# Test linear SVM
test_that("Simple svm predicts correctly", {
  model <- svm(X, y, C=1)
  y_pred <- model(X)
  expect_equal(y, y_pred)
})

test_that("Linear kernel svm predicts correctly", {
  ckernel <- tuned_kernel(lin_kernel)
  model <- kernel_svm(X, y, C=1, ckernel)
  y_pred <- model(X)
  expect_equal(y, y_pred)
})

# Create point clouds with non linear boundary
n <- 50
c1_x1 <- rnorm(n, 6, 1)
c1_x2 <- rnorm(n, 6, 1)

c2_x1 <- c(runif(n/2, 1, 10), runif(n/2, 1, 2))
c2_x2 <- c(runif(n/2, 1, 2), runif(n/2, 1, 10))

X <- matrix(c(c1_x1, c2_x1, c1_x2, c2_x2), ncol=2)
y <- c(rep(-1, n), rep(1, n))

# Test kernel SVM
test_that("Polynomial kernel svm predicts > 95% correctly (data scaled first)", {
  Xs <- scale_dat(X, X)[,2:3]
  ckernel <- tuned_kernel(poly_kernel, b=3)
  model <- kernel_svm(Xs, y, C=1, ckernel)
  y_pred <- model(Xs)
  expect_true(sum(y == y_pred)/length(y) > 0.95)
})

test_that("RBF kernel svm predicts correctly > 95% correctly", {
  ckernel <- tuned_kernel(rbf_kernel, sigma=5)
  model <- kernel_svm(X, y, C=1, ckernel)
  y_pred <- model(X)
  expect_true(sum(y == y_pred)/length(y) > 0.95)
})

test_that("Trigonometric kernel svm predicts correctly > 90% correctly (more variable)", {
  ckernel <- tuned_kernel(trig_kernel, b=5)
  model <- kernel_svm(X, y, C=1, ckernel)
  y_pred <- model(X)
  expect_true(sum(y == y_pred)/length(y) > 0.9)
})

# plot coloured by true and predicted class
# plot(X[,1], X[,2], col=as.factor(y))
# points(X[,1], X[,2], col=as.factor(y_pred), pch=16)

