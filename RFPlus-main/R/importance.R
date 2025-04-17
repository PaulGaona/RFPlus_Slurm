#' Compute permutation importance for RF+ estimators.
#'
#' @inheritParams rfplus
#' @param object An object of class \code{rfplus}.
#' @param metric A function that takes two vectors and returns a scalar.
#' @param grouped_features A list of vectors of feature names. Each vector
#'   corresponds to a group of features that should be permuted together.
#' @param B The number of permutations to perform. If \code{B} is a list of
#'   vectors of indices, then the indices will be used to permute the features.
#' @param local Whether to compute local feature importance.
#' @param local_metric A function that takes two vectors and returns a vector
#'   of the same length.
#' @param return_preds Whether to return the predictions for each permutation.
#'
#' @export
permute_fi <- function(object, x, y,
                       metric = yardstick::rmse_vec,
                       grouped_features = NULL,
                       B = 10,
                       local = FALSE,
                       local_metric = metric,
                       return_preds = FALSE) {

  if (!is.list(B)) {
    permute_idxs <- purrr::map(1:B, ~ sample(1:nrow(x)))
  } else {
    permute_idxs <- B
  }

  if (is.null(grouped_features)) {
    grouped_features <- colnames(x)
    names(grouped_features) <- colnames(x)
  }

  orig_preds <- predict_tree(object, x = x)
  orig_score <- metric(truth = y, estimate = orig_preds)
  if (local) {
    local_orig_score <- purrr::map2_dbl(y, orig_preds, local_metric)
  }

  fi_preds <- purrr::map(
    grouped_features,
    function(j) {
      if (is.null(j)) {
        return(purrr::map(permute_idxs, ~ orig_preds))
      }
      partial_preds <- purrr::map(
        permute_idxs,
        function(permute_idx) {
          x_mod <- x
          x_mod[, j] <- x_mod[permute_idx, j]
          predict_tree(object, x = x_mod)
        }
      )
    }
  )
  fis <- purrr::map(
    fi_preds,
    function(partial_preds) {
      score <- purrr::map_dbl(
        partial_preds, ~ metric(truth = y, estimate = .x) - orig_score
      ) %>%
        mean()
    }
  ) %>%
    tibble::as_tibble()

  if (local) {
    local_fis <- purrr::map(
      fi_preds,
      function(partial_preds) {
        score <- purrr::map(
          partial_preds,
          ~ purrr::map2_dbl(y, .x, local_metric) - local_orig_score
        ) %>%
          purrr::reduce(`+`) / length(partial_preds)
      }
    ) %>%
      tibble::as_tibble()
  }

  if (return_preds) {
    if (local) {
      return(list(global = fis, local = local_fis, preds = fi_preds))
    } else {
      return(list(global = fis, preds = fi_preds))
    }
  } else {
    if (local) {
      return(list(global = fis, local = local_fis))
    } else {
      return(fis)
    }
  }
}


#' Compute MDI+ importance for RF+ estimators.
#'
#' @inheritParams rfplus
#' @inheritParams permute_fi
#' @param col_means A vector of column means to use for the global feature
#'   importance. If \code{NULL}, then the column means of \code{x} will be used.
#'
#' @export
mdiplus_fi <- function(object, x, y, col_means = NULL,
                       metric = rsq_vec_narm, grouped_features = NULL,
                       local = FALSE, local_metric = yardstick::rmse_vec) {
  if (is.null(col_means)) {
    col_means <- apply(x, 2, mean)
  }

  if (is.null(grouped_features)) {
    grouped_features <- colnames(x)
    names(grouped_features) <- colnames(x)
  }

  fis <- purrr::map(
    grouped_features,
    function(j) {
      if (is.null(j)) {
        partial_preds <- predict_tree(
          object, x = matrix(col_means, nrow = 1)
        ) %>%
          rep(nrow(x))
      } else {
        mean_nodes <- setdiff(colnames(x), j)
        x_mod <- x
        if (length(mean_nodes) == 1) {
          x_mod[, mean_nodes] <- col_means[mean_nodes]
        } else {
          x_mod[, mean_nodes] <- matrix(
            col_means[mean_nodes],
            nrow = nrow(x), ncol = length(mean_nodes), byrow = TRUE
          )
        }
        partial_preds <- predict_tree(object, x = x_mod)
      }
      score <- metric(truth = y, estimate = partial_preds)
      if (local) {
        local_score <- purrr::map2_dbl(y, partial_preds, local_metric)
      } else {
        local_score <- NULL
      }
      return(list(global = score, local = local_score))
    }
  )

  global_fis <- purrr::map(fis, "global") %>%
    tibble::as_tibble()
  if (local) {
    local_fis <- purrr::map(fis, "local") %>%
      tibble::as_tibble()
  }

  if (local) {
    return(list(global = global_fis, local = local_fis))
  } else {
    return(global_fis)
  }
}


#' Compute feature importance for RF+ estimators.
#'
#' @inheritParams rfplus
#' @inheritParams permute_fi
#' @param type The type of feature importance to compute. Can be one of
#'   \code{"permute"} or \code{"mdi+"}. See \code{\link{permute_fi}} and
#'   \code{\link{mdiplus_fi}} for details.
#'
#' @export
get_feature_importances <- function(object, x, y, type = c("permute", "mdi+"),
                                    B = 10, local = FALSE, use_rcpp = TRUE) {
  type <- match.arg(type)

  if (is.matrix(x)) {
    x <- data.frame(x)
  }

  rf_fit <- object$rf_fit
  rfplus_fits <- object$rfplus_fits
  tree_infos <- object$tree_infos
  # forest_paths <- fit$forest_paths
  x_train <- object$x_train
  # psis_train <- fit$psis_train
  x_scale_factors <- object$x_scale_factors
  include_raw <- object$include_raw
  normalize_stumps <- object$normalize_stumps
  normalize_raw <- object$normalize_raw
  sample_split <- object$sample_split

  if (identical(type, "permute")) {
    # TODO: use different metric for classification
    metric <- yardstick::rmse_vec
  }

  train_test_psis <- get_train_test_psis(object, x = x)
  psis_train <- train_test_psis$train
  psis <- train_test_psis$test

  permute_idxs <- purrr::map(1:B, ~ sample(1:nrow(x)))

  fi_scores_ls <- purrr::pmap(
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

      grouped_features <- get_grouped_tree_features(
        colnames(x), colnames(x_aug), tree_info
      )

      # get feature importance
      if (identical(type, "mdi+")) {
        fis <- mdiplus_fi(
          object = rfplus_fit, x = x_aug, y = y,
          col_means = apply(x_train_aug, 2, mean),
          metric = yardstick::rsq_vec,
          grouped_features = grouped_features,
          local = local, local_metric = yardstick::rmse_vec
        )
      } else if (identical(type, "permute")) {
        fis <- permute_fi(
          object = rfplus_fit, x = x_aug, y = y,
          metric = yardstick::rmse_vec,
          grouped_features = grouped_features,
          B = permute_idxs, local = FALSE, return_preds = TRUE
        )$preds
      }
    }
  )

  if (identical(type, "permute")) {
    metric <- yardstick::rmse_vec
    local_metric <- metric
    orig_preds <- predict(object, x)
    orig_score <- metric(truth = y, estimate = orig_preds)
    fi_scores <- purrr::map(
      purrr::list_transpose(fi_scores_ls),
      function(var_ls) {
        score <- purrr::map_dbl(
          purrr::list_transpose(var_ls),
          function(.x) {
            preds <- purrr::reduce(.x, `+`) / length(.x)
            return(metric(truth = y, estimate = preds) - orig_score)
          }
        ) %>%
          mean()
      }
    ) %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "var", values_to = "importance"
      )
    if (local) {
      local_orig_score <- purrr::map2_dbl(y, orig_preds, local_metric)
      local_fi_scores <- purrr::map(
        purrr::list_transpose(fi_scores_ls),
        function(var_ls) {
          score <- purrr::map(
            purrr::list_transpose(var_ls),
            function(.x) {
              preds <- purrr::reduce(.x, `+`) / length(.x)
              return(purrr::map2_dbl(y, preds, local_metric) - local_orig_score)
            }
          ) %>%
            purrr::reduce(`+`) / length(purrr::list_transpose(var_ls))
        }
      ) %>%
        tibble::as_tibble()
    }
  } else if (identical(type, "mdi+")) {
    if (local) {
      global_fi_scores_ls <- purrr::map(fi_scores_ls, "global")
    } else {
      global_fi_scores_ls <- fi_scores_ls
    }
    fi_scores <- global_fi_scores_ls %>%
      purrr::list_rbind() %>%
      dplyr::summarise(
        dplyr::across(
          tidyselect::everything(), ~ mean(.x, na.rm = TRUE)
        )
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "var", values_to = "importance"
      )
    if (local) {
      local_fi_scores <- purrr::map(fi_scores_ls, "local") %>%
        purrr::reduce(`+`) / length(fi_scores_ls)
    }
  }

  if (local) {
    return(list(global = fi_scores, local = local_fi_scores))
  } else {
    return(fi_scores)
  }
}


#' Get stump features corresponding to each raw feature
#'
#' @param orig_colnames Original (raw) column names.
#' @param aug_colnames Augmented (raw + stump) column names.
#' @param tree_info Output of [ranger::treeInfo()] for a single tree.
#'
#' @export
get_grouped_tree_features <- function(orig_colnames, aug_colnames, tree_info) {
  grouped_features <- purrr::map(
    orig_colnames,
    function(j) {
      j_nodes <- tree_info %>%
        dplyr::filter(splitvarName == !!j) %>%
        dplyr::pull(nodeID)
      if (length(j_nodes) == 0) {
        return(NULL)
      } else {
        j_nodes <- paste0(".node", j_nodes)
      }
      xj_features <- aug_colnames[
        stringr::str_detect(aug_colnames, paste0("^", j, "$")) |
          stringr::str_detect(aug_colnames, paste0("^", j, "\\.[0-9]*\\.\\."))
      ]
      c(xj_features, j_nodes)
    }
  ) %>%
    setNames(orig_colnames)
  return(grouped_features)
}
