test_that("rfplus works", {

  n <- 100
  n_test <- 50
  p <- 10
  s <- 2
  noise_sd <- 0.1
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- matrix(1, nrow = s, ncol = 1)
  y <- c(x[, 1:s, drop = FALSE] %*% beta + rnorm(n, sd = noise_sd))
  x_test <- matrix(rnorm(n_test * p), nrow = n, ncol = p)
  y_test <- c(x_test[, 1:s, drop = FALSE] %*% beta + rnorm(n, sd = noise_sd))

  rfplus_fit <- rfplus(x = x, y = y)
  rfplus_preds <- predict(rfplus_fit, x = x_test)

  rfplus_ridge_fit <- rfplus(x = x, y = y, lambda_x = 1, lambda_t = 2)
  rfplus_ridge_preds <- predict(rfplus_ridge_fit, x = x_test)

  lambda_xs <- c(1, 10)
  lambda_ts <- c(2, 4)
  rfplus_cv_fit <- rfplus_cv(
    x = x, y = y, lambda_xs = lambda_xs, lambda_ts = lambda_ts
  )
  rfplus_cv_preds <- predict(rfplus_cv_fit, x = x_test)

  mdiplus_fis <- get_feature_importances(
    rfplus_cv_fit, x = x_test, y = y_test, type = "mdi+"
  )
  # mdiplus_fis
  # expect_true(
  #   all(rank(-mdiplus_fis$importance)[1:s] <= s)
  # )
  permute_fis <- get_feature_importances(
    rfplus_cv_fit, x = x_test, y = y_test, type = "permute"
  )
  # permute_fis
  # expect_true(
  #   all(rank(permute_fis$importance)[1:s] <= s)
  # )

  mdiplus_local_fis <- get_feature_importances(
    rfplus_cv_fit, x = x_test, y = y_test, type = "mdi+", local = TRUE
  )
  permute_local_fis <- get_feature_importances(
    rfplus_cv_fit, x = x_test, y = y_test, type = "permute", local = TRUE
  )
})
