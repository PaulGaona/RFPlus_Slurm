#' Random Forest Plus (RF+)
#'
#' @param x A matrix or data frame of predictors.
#' @param y A vector of responses.
#' @param include_raw Whether to include the raw features in the RF+ estimator.
#' @param normalize_stumps Whether to normalize the stumps in the RF+ estimator.
#' @param normalize_raw Whether to normalize the raw features in the RF+
#'   estimator.
#' @param family The family of the response. Either "linear" or "logistic".
#' @param sample_split How to perform sample splitting for the RF+ estimator.
#'   Must be one of "none" or "oob".
#' @param ntrees The number of trees to grow.
#' @param mtry The number of variables to sample as candidates at each split.
#' @param num.threads The number of threads to use for training RF+.
#' @param lambda_x The regularization parameter for the raw features.
#' @param lambda_t The regularization parameter for the stumps.
#' @param parallel Whether to use parallel processing.
#' @param low_mem Whether to use low memory mode. Should remove this argument
#'   before release. Keep only for debugging.
#' @param use_rcpp Whether to use Rcpp. Should remove this argument before
#'   release. Keep only for debugging.
#' @param ... Other arguments passed to \code{ranger::ranger}.
#'
#' @export
rfplus <- function(x, y, include_raw = TRUE,
                   normalize_stumps = FALSE, normalize_raw = FALSE,
                   family = c("linear", "logistic"),
                   sample_split = c("none", "oob"),
                   ntrees = 500, mtry = NULL, num.threads = 1,
                   lambda_x = 0, lambda_t = lambda_x,
                   parallel = FALSE,
                   # get rid of arguments later
                   low_mem = TRUE, use_rcpp = TRUE, ...) {
  family <- match.arg(family)
  sample_split <- match.arg(sample_split)
  if (is.null(mtry)) {
    if (family == "linear") {
      mtry <- ncol(x) / 3
    } else if (family == "logistic") {
      mtry <- sqrt(ncol(x))
    }
  }
  if (is.matrix(x)) {
    x <- data.frame(x)
  }
  if (family == "linear") {
    classification <- FALSE
  } else if (family == "logistic") {
    classification <- TRUE
  }
  if (parallel) {
    map_fun <- furrr::future_map
    imap_fun <- furrr::future_imap
    pmap_fun <- furrr::future_pmap
  } else {
    map_fun <- purrr::map
    imap_fun <- purrr::imap
    pmap_fun <- purrr::pmap
  }
  if (length(lambda_x) == 1) {
    lambda_x <- rep(lambda_x, ntrees)
  } else if (length(lambda_x) != ntrees) {
    stop("lambda_x must be a scalar or vector of length ntrees.")
  }
  if (length(lambda_t) == 1) {
    lambda_t <- rep(lambda_t, ntrees)
  } else if (length(lambda_t) != ntrees) {
    stop("lambda_t must be a scalar or vector of length ntrees.")
  }

  train_df <- dplyr::bind_cols(x, .y = y)
  rf_fit <- ranger::ranger(
    formula = .y ~ .,
    data = train_df,
    num.trees = ntrees,
    mtry = mtry,
    num.threads = num.threads,
    keep.inbag = TRUE,
    classification = classification,
    node.stats = TRUE,
    ...
  )

  tree_infos <- purrr::map(1:ntrees, ~ ranger::treeInfo(rf_fit, .x))
  node_preds <- predict(
    rf_fit, x, type = "terminalNodes", num.threads = num.threads
  )$predictions
  if (use_rcpp) {
    forest_paths <- get_forest_paths_fast(tree_infos)
  } else {
    forest_paths <- get_forest_paths(tree_infos)
  }

  unordered_factors <- colnames(x)[!rf_fit$forest$is.ordered]
  x_psi <- get_x_psi(x, x, rf_fit)
  psis <- imap_fun(
    tree_infos,
    function(tree_info, tree_idx) {
      if (use_rcpp) {
        psi <- extract_psi_fast(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      } else {
        psi <- extract_psi(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      }
    }
  )

  rfplus_fits <- pmap_fun(
    list(
      psi = psis,
      tree_info = tree_infos,
      inbag_counts = rf_fit$inbag.counts,
      lam_x = lambda_x,
      lam_t = lambda_t
    ),
    function(psi, tree_info, inbag_counts, lam_x, lam_t) {
      x_train <- get_augmented_x(
        x = x, psi = psi, tree_info = tree_info,
        include_raw = include_raw, normalize = normalize_raw
      )
      y_train <- y
      if (sample_split == "oob") {
        x_train <- x_train[inbag_counts > 0, , drop = FALSE]
        y_train <- y_train[inbag_counts > 0]
      }
      if (include_raw) {
        ridge_lambda <- c(
          rep(lam_x, length.out = ncol(x_train) - ncol(psi)),
          rep(lam_t, length.out = ncol(psi))
        )
      } else {
        ridge_lambda <- rep(lam_t, length.out = ncol(x_train))
      }
      if (all(lam_x == 0)) {
        aug_train_df <- dplyr::bind_cols(x_train, .y = y_train)
        if (classification) {
          fit <- glm(.y ~ ., data = aug_train_df, family = "binomial")
        } else {
          fit <- lm(.y ~ ., data = aug_train_df)
        }
      } else {
        if (classification) {
          family <- "binomial"
        } else {
          family <- "gaussian"
        }
        fit <- glmnet::glmnet(
          x = x_train, y = y_train, alpha = 0,
          lambda = sum(ridge_lambda) / ncol(x_train),
          penalty.factor = ridge_lambda,
          family = family
        )
      }
      attr(fit, "x_scale_factors") <- attr(x_train, "x_scale_factors")
      return(fit)
    }
  )

  x_scale_factors <- purrr::map(rfplus_fits, ~ attr(.x, "x_scale_factors"))

  out <- list(
    rf_fit = rf_fit,
    rfplus_fits = rfplus_fits,
    tree_infos = tree_infos,
    forest_paths = forest_paths,
    x_train = x,
    psis_train = psis,
    x_scale_factors = x_scale_factors,
    lambda_x = lambda_x,
    lambda_t = lambda_t,
    include_raw = include_raw,
    normalize_stumps = normalize_stumps,
    normalize_raw = normalize_raw,
    classification = classification,
    sample_split = sample_split,
    node_preds = node_preds
  )
  class(out) <- "rfplus"

  return(out)
}

#' Random Forest Plus Dummy (RF+)
#'
#' @inheritParams rfplus
#'

#' @export
rfplus.dummy <- function(x, y, include_raw = TRUE,
                   normalize_stumps = FALSE, normalize_raw = FALSE,
                   family = c("linear", "logistic"),
                   sample_split = c("none", "oob"),
                   ntrees = 500, mtry = NULL, num.threads = 1,
                   lambda_x = 0, lambda_t = lambda_x,
                   parallel = FALSE,
                   # get rid of arguments later
                   low_mem = TRUE, use_rcpp = TRUE, ...) {
  family <- match.arg(family)
  sample_split <- match.arg(sample_split)
  if (is.null(mtry)) {
    if (family == "linear") {
      mtry <- ncol(x) / 3
    } else if (family == "logistic") {
      mtry <- sqrt(ncol(x))
    }
  }
  if (is.matrix(x)) {
    x <- data.frame(x)
  }
  if (family == "linear") {
    classification <- FALSE
  } else if (family == "logistic") {
    classification <- TRUE
  }
  if (parallel) {
    map_fun <- furrr::future_map
    imap_fun <- furrr::future_imap
    pmap_fun <- furrr::future_pmap
  } else {
    map_fun <- purrr::map
    imap_fun <- purrr::imap
    pmap_fun <- purrr::pmap
  }
  if (length(lambda_x) == 1) {
    lambda_x <- rep(lambda_x, ntrees)
  } else if (length(lambda_x) != ntrees) {
    stop("lambda_x must be a scalar or vector of length ntrees.")
  }
  if (length(lambda_t) == 1) {
    lambda_t <- rep(lambda_t, ntrees)
  } else if (length(lambda_t) != ntrees) {
    stop("lambda_t must be a scalar or vector of length ntrees.")
  }
  
  train_df <- dplyr::bind_cols(x, .y = y)
  rf_fit <- ranger::ranger(
    formula = .y ~ .,
    data = train_df,
    num.trees = ntrees,
    mtry = mtry,
    num.threads = num.threads,
    keep.inbag = TRUE,
    classification = classification,
    ...
  )
  
  tree_infos <- purrr::map(1:ntrees, ~ ranger::treeInfo(rf_fit, .x))
  node_preds <- predict(
    rf_fit, x, type = "terminalNodes", num.threads = num.threads
  )$predictions
  if (use_rcpp) {
    forest_paths <- get_forest_paths_fast(tree_infos)
  } else {
    forest_paths <- get_forest_paths(tree_infos)
  }
  
  unordered_factors <- colnames(x)[!rf_fit$forest$is.ordered]
  x_psi <- get_x_psi(x, x, rf_fit)
  psis <- imap_fun(
    tree_infos,
    function(tree_info, tree_idx) {
      if (use_rcpp) {
        psi <- extract_psi_fast(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      } else {
        psi <- extract_psi(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      }
    }
  )
  
  rfplus_fits <- pmap_fun(
    list(
      psi = psis,
      tree_info = tree_infos,
      inbag_counts = rf_fit$inbag.counts,
      lam_x = lambda_x,
      lam_t = lambda_t
    ),
    function(psi, tree_info, inbag_counts, lam_x, lam_t) {
      x_train <- get_augmented_x(
        x = x, psi = psi, tree_info = tree_info,
        include_raw = include_raw, normalize = normalize_raw
      )
      y_train <- y
      if (sample_split == "oob") {
        x_train <- x_train[inbag_counts > 0, , drop = FALSE]
        y_train <- y_train[inbag_counts > 0]
      }
      if (include_raw) {
        ridge_lambda <- c(
          rep(lam_x, length.out = ncol(x_train) - ncol(psi)),
          rep(lam_t, length.out = ncol(psi))
        )
      } else {
        ridge_lambda <- rep(lam_t, length.out = ncol(x_train))
      }
      if (all(lam_x == 0)) {
        aug_train_df <- dplyr::bind_cols(x_train, .y = y_train)
        if (classification) {
          fit <- glm(.y ~ ., data = aug_train_df, family = "binomial")
        } else {
          fit <- lm(.y ~ ., data = aug_train_df)
        }
      } else {
        if (classification) {
          family <- "binomial"
        } else {
          family <- "gaussian"
        }
        fit <- glmnet::glmnet(
          x = x_train, y = y_train, alpha = 0,
          lambda = sum(ridge_lambda) / ncol(x_train),
          penalty.factor = ridge_lambda,
          family = family
        )
      }
      attr(fit, "x_scale_factors") <- attr(x_train, "x_scale_factors")
      return(fit)
    }
  )
  
  x_scale_factors <- purrr::map(rfplus_fits, ~ attr(.x, "x_scale_factors"))
  
  out <- list(
    rf_fit = rf_fit,
    rfplus_fits = rfplus_fits,
    tree_infos = tree_infos,
    forest_paths = forest_paths,
    x_train = x,
    psis_train = psis,
    x_scale_factors = x_scale_factors,
    lambda_x = lambda_x,
    lambda_t = lambda_t,
    include_raw = include_raw,
    normalize_stumps = normalize_stumps,
    normalize_raw = normalize_raw,
    classification = classification,
    sample_split = sample_split,
    node_preds = node_preds
  )
  class(out) <- "rfplus"
  
  return(out)
}

#' Random Forest Plus (RF+) with cross-validation wrapepr.
#'
#' @inheritParams rfplus
#' @param ntrees_cv Number of trees to use for cross-validation. Defaults to
#'   \code{ntrees}. If \code{ntrees_cv} is less than \code{ntrees}, then
#'   \code{ntrees_cv} trees will be tuned via cross-validation and the remaining
#'   \code{ntrees - ntrees_cv} trees will choose a hyperparameter set randomly
#'   from the best hyperparameters from the cross-validated \code{ntrees_cv}
#'   trees.
#' @param lambda_xs Vector or list of regularization parameters for the raw
#'   features to search over.
#' @param lambda_ts Vector or list of regularization parameters for the stump
#'   features to search over.
#' @param cv Number of folds for cross-validation.
#'
#' @export
rfplus_cv <- function(x, y, include_raw = TRUE,
                      normalize_stumps = FALSE, normalize_raw = FALSE,
                      family = c("linear", "logistic"),
                      sample_split = c("none", "oob"),
                      ntrees = 500, ntrees_cv = ntrees,
                      mtry = NULL, num.threads = 1,
                      lambda_xs = 0, lambda_ts = NULL,
                      cv = 5, parallel = FALSE,
                      # get rid of arguments later
                      low_mem = TRUE, use_rcpp = TRUE, ...) {
  family <- match.arg(family)
  sample_split <- match.arg(sample_split)
  if (is.null(mtry)) {
    if (family == "linear") {
      mtry <- ncol(x) / 3
    } else if (family == "logistic") {
      mtry <- sqrt(ncol(x))
    }
  }
  if (is.matrix(x)) {
    x <- data.frame(x)
  }
  if (family == "linear") {
    classification <- FALSE
  } else if (family == "logistic") {
    classification <- TRUE
  }
  if (parallel) {
    map_fun <- furrr::future_map
    imap_fun <- furrr::future_imap
    pmap_fun <- furrr::future_pmap
  } else {
    map_fun <- purrr::map
    imap_fun <- purrr::imap
    pmap_fun <- purrr::pmap
  }

  train_df <- dplyr::bind_cols(x, .y = y)
  rf_fit <- ranger::ranger(
    formula = .y ~ .,
    data = train_df,
    num.trees = ntrees,
    mtry = mtry,
    num.threads = num.threads,
    keep.inbag = TRUE,
    classification = classification,
    ...
  )

  tree_infos <- purrr::map(1:ntrees, ~ ranger::treeInfo(rf_fit, .x))
  node_preds <- predict(
    rf_fit, x, type = "terminalNodes", num.threads = num.threads
  )$predictions
  if (use_rcpp) {
    forest_paths <- get_forest_paths_fast(tree_infos)
  } else {
    forest_paths <- get_forest_paths(tree_infos)
  }

  unordered_factors <- colnames(x)[!rf_fit$forest$is.ordered]
  x_psi <- get_x_psi(x, x, rf_fit)
  psis <- imap_fun(
    tree_infos,
    function(tree_info, tree_idx) {
      if (use_rcpp) {
        psi <- extract_psi_fast(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      } else {
        psi <- extract_psi(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      }
    }
  )


  if (is.null(lambda_ts)) {
    param_grid <- tibble::tibble(
      lambda_x = lambda_xs,
      lambda_t = lambda_xs
    )
  } else {
    param_grid <- tidyr::expand_grid(
      lambda_x = lambda_xs,
      lambda_t = lambda_ts
    )
  }

  cv_losses <- NULL
  best_cv_params <- NULL
  rfplus_fits <- map_fun(
    1:ntrees,
    function(i) {
      psi <- psis[[i]]
      tree_info <- tree_infos[[i]]
      inbag_counts <- rf_fit$inbag.counts[[i]]
      x_train <- get_augmented_x(
        x = x, psi = psi, tree_info = tree_info,
        include_raw = include_raw, normalize = normalize_raw
      )
      y_train <- y
      if (sample_split == "oob") {
        x_train <- x_train[inbag_counts > 0, , drop = FALSE]
        y_train <- y_train[inbag_counts > 0]
      }
      if (all(purrr::map_lgl(lambda_xs, ~ all(.x == 0)))) {
        aug_train_df <- dplyr::bind_cols(x_train, .y = y_train)
        if (classification) {
          fit <- glm(.y ~ ., data = aug_train_df, family = "binomial")
        } else {
          fit <- lm(.y ~ ., data = aug_train_df)
        }
      } else {
        if (classification) {
          family <- "binomial"
          response_type <- "class"
          metric_fun <- yardstick::accuracy_vec
        } else {
          family <- "gaussian"
          response_type <- "response"
          metric_fun <- function(...) -yardstick::rmse_vec(...)
        }
        if (is.list(lambda_xs)) {
          cv_losses <- matrix(NA, nrow = nrow(param_grid), ncol = cv)
          cv_foldids <- sample(1:cv, size = nrow(x_train), replace = TRUE)
          for (k in 1:cv) {
            train_idx <- cv_foldids != k
            valid_idx <- cv_foldids == k
            # TODO: speed up by using glmnet::cv.glmnet(?) or manual eigendecomposition and use ntrees_cv
            for (param_idx in 1:nrow(param_grid)) {
              cur_lambdas <- c(
                rep(param_grid$lambda_x[[param_idx]],
                    length.out = ncol(x_train) - ncol(psi)),
                rep(param_grid$lambda_t[[param_idx]],
                    length.out = ncol(psi))
              )
              fold_fit <- glmnet::glmnet(
                x = x_train[train_idx, , drop = FALSE],
                y = y_train[train_idx],
                alpha = 0,
                # to fix glmnet internal rescaling of penalty factor
                lambda = sum(cur_lambdas) / ncol(x_train),
                penalty.factor = cur_lambdas,
                family = family
              )
              fold_preds <- predict(
                fold_fit, x_train[valid_idx, , drop = FALSE],
                type = response_type
              )[, 1]
              cv_losses[param_idx, k] <- metric_fun(
                y_train[valid_idx], fold_preds
              )
            }
          }
          best_cv_params <- param_grid[which.max(rowMeans(cv_losses)), ]
          best_lambdas <- c(
            rep(best_cv_params$lambda_x[[1]],
                length.out = ncol(x_train) - ncol(psi)),
            rep(best_cv_params$lambda_t[[1]],
                length.out = ncol(psi))
          )
          fit <- glmnet::glmnet(
            x = x_train, y = y_train,
            alpha = 0,
            lambda = sum(best_lambdas) / ncol(x_train),
            penalty.factor = best_lambdas,
            family = family
          )
        } else {
          # TODO: fix/add CV for lambda_t
          cv_fit <- glmnet::cv.glmnet(
            x = x_train, y = y_train,
            alpha = 0, lambda = lambda_xs, family = family
          )
          fit <- glmnet::glmnet(
            x = x_train, y = y_train, alpha = 0, lambda = cv_fit$lambda.1se,
            family = family
          )
        }
      }

      attr(fit, "x_scale_factors") <- attr(x_train, "x_scale_factors")
      return(fit)
    }
  )

  x_scale_factors <- purrr::map(rfplus_fits, ~ attr(.x, "x_scale_factors"))

  out <- list(
    rf_fit = rf_fit,
    rfplus_fits = rfplus_fits,
    cv_losses = cv_losses,
    best_cv_params = best_cv_params,
    tree_infos = tree_infos,
    forest_paths = forest_paths,
    x_train = x,
    psis_train = psis,
    x_scale_factors = x_scale_factors,
    ntrees_cv = ntrees_cv,
    include_raw = include_raw,
    normalize_stumps = normalize_stumps,
    normalize_raw = normalize_raw,
    classification = classification,
    sample_split = sample_split
  )
  class(out) <- c("rfplus", "rfplus_cv")

  return(out)
}


#' @keywords internal
predict_tree <- function(object, x, type = c("response")) {
  type <- match.arg(type)
  if ("glmnet" %in% class(object)) {
    rfplus_preds <- predict(object, x, type = type)
  } else {
    rfplus_preds <- predict(object, data.frame(x), type = type)
  }
  if (isTRUE(ncol(rfplus_preds) == 1)) {
    rfplus_preds <- c(rfplus_preds)
  }
  return(rfplus_preds)
}


#' Predict using Random Forest Plus (RF+)
#'
#' @param object An object of class \code{rfplus}.
#' @param x A matrix or data frame of observations.
#' @param type The type of prediction to return.
#' @param return_all Whether to return predictions from all trees (rather than
#'   the average prediction).
#' @param use_rcpp Whether to use the Rcpp. Should remove this argument before
#'   release. Keep only for debugging.
#'
#' @export
predict.rfplus <- function(object, x, type = c("response"),
                           return_all = FALSE,
                           # get rid of arguments later
                           use_rcpp = TRUE) {
  type <- match.arg(type)
  if (is.matrix(x)) {
    x <- data.frame(x)
  }

  rf_fit <- object$rf_fit
  rfplus_fits <- object$rfplus_fits
  tree_infos <- object$tree_infos
  forest_paths <- object$forest_paths
  x_train <- object$x_train
  psis_train <- object$psis_train
  x_scale_factors <- object$x_scale_factors
  include_raw <- object$include_raw
  normalize_stumps <- object$normalize_stumps
  normalize_raw <- object$normalize_raw
  sample_split <- object$sample_split

  train_test_psis <- get_train_test_psis(object, x = x, use_rcpp = use_rcpp)
  psis_train <- train_test_psis$train
  psis <- train_test_psis$test

  preds <- purrr::pmap(
    list(
      psi = psis,
      psi_train = psis_train,
      tree_info = tree_infos,
      x_scale_factor = x_scale_factors,
      rfplus_fit = rfplus_fits,
      inbag_counts = rf_fit$inbag.counts
    ),
    function(psi, psi_train, tree_info, x_scale_factor, rfplus_fit, inbag_counts) {
      x_train_aug <- get_augmented_x(
        x = x_train, psi = psi_train, tree_info = tree_info,
        include_raw = include_raw, normalize = normalize_raw,
        normalize_factor = x_scale_factor
      )
      if (sample_split == "oob") {
        x_train_aug <- x_train_aug[inbag_counts > 0, , drop = FALSE]
      }
      x_aug <- get_augmented_x(
        x = x, psi = psi, tree_info = tree_info,
        include_raw = include_raw, normalize = normalize_raw,
        normalize_factor = x_scale_factor
      )
      x_aug_full <- rbind(x_train_aug, x_aug)
      rfplus_preds <- predict_tree(rfplus_fit, x = x_aug, type = type)
    }
  )
  if (!return_all) {
    preds <- purrr::reduce(preds, `+`) / rf_fit$num.trees
    if (isTRUE(ncol(preds) == 1)) {
      preds <- c(preds)
    }
  }
  return(preds)
}
