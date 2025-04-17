# Note: Many of these functions should not be exported in the release version.
# Exporting now for debugging purposes.


#' Extract all root-to-leaf paths in a forest.
#'
#' @param tree_infos List of size num.trees with each entry being the output of
#'   [ranger::treeInfo()].
#'
#' @returns A list of size num.trees with each entry being a list of decision
#'   paths in each tree.
#' @export
get_forest_paths_fast <- function(tree_infos) {
  purrr::map(tree_infos, ~get_tree_paths_fast(.x))
}


#' Extract all root-to-leaf paths in a forest (non-Rcpp version).
#'
#' @inheritParams get_forest_paths_fast
#' @export
get_forest_paths <- function(tree_infos) {
  purrr::map(tree_infos, ~get_tree_paths(.x))
}

#' Extract all root-to-leaf paths in a tree.
#'
#' @param tree_info Output of [ranger::treeInfo()] for a single tree.
#'
#' @returns A list of decision paths in a single tree.
#' @export
get_tree_paths_fast <- function(tree_info) {
  terminal_node_ids <- tree_info$nodeID[tree_info$terminal]
  inner_tree_info <- tree_info %>%
    dplyr::filter(!terminal)
  tree_paths <- get_tree_paths_cpp(
    terminal_node_ids = terminal_node_ids,
    left_child_ids = inner_tree_info$leftChild,
    right_child_ids = inner_tree_info$rightChild,
    node_ids = inner_tree_info$nodeID
  ) %>%
    setNames(as.character(terminal_node_ids))
  return(tree_paths)
}


#' Extract all root-to-leaf paths in a tree (non-Rcpp version).
#'
#' @inheritParams get_tree_paths_fast
#' @export
get_tree_paths <- function(tree_info) {

  terminal_node_ids <- tree_info$nodeID[tree_info$terminal]
  inner_tree_info <- tree_info %>%
    dplyr::filter(!terminal)
  tree_paths <- list()
  for (terminal_node_id in terminal_node_ids) {
    node_id <- terminal_node_id
    tree_path <- c()
    while (node_id != 0) {
      if (node_id %in% inner_tree_info$leftChild) {
        idx <- inner_tree_info$leftChild == node_id
      } else {
        idx <- inner_tree_info$rightChild == node_id
      }
      node_id <- inner_tree_info$nodeID[idx]
      tree_path <- c(tree_path, node_id)
    }
    tree_paths[[as.character(terminal_node_id)]] <- tree_path
  }
  return(tree_paths)
}


#' Extract Psi (i.e., stump) matrix from tree.
#'
#' @inheritParams get_tree_paths_fast
#' @param x Data matrix for which to compute Psi.
#' @param tree_paths Output of [get_tree_paths_fast()].
#' @param node_preds Terminal node predictions corresponding to x.
#' @param psi_train (Optional) Psi matrix for training data. If provided, the
#'   output will be normalized using the same normalizing constants as the
#'   training data. Otherwise, the output will be normalized using the
#'   normalizing constants computed from the input data. Only used if
#'   \code{normalize = TRUE}.
#' @param inbag_counts (Optional) Number of times each training observation was
#'   sampled in the tree (from RF+). Only used if \code{normalize = TRUE}.
#' @param normalize Whether to normalize the stumps so that larger stumps
#'   have larger magnitudes.
#' @param unordered_factors Names of unordered factors in x.
#'
#' @export
extract_psi_fast <- function(x, tree_info, tree_paths, node_preds,
                             psi_train = NULL, inbag_counts = NULL,
                             normalize = FALSE, unordered_factors = NULL) {
  x_mat <- as.matrix(x)
  x_colnames_df <- tibble::tibble(
    name = colnames(x),
    idx = 0:(length(colnames(x)) - 1)
  )

  inner_tree_info <- tree_info %>%
    dplyr::filter(!terminal) %>%
    dplyr::left_join(x_colnames_df, by = c("splitvarName" = "name"))

  if (length(unordered_factors) == 0) {
    psi <- extract_psi_cpp(
      x = x_mat,
      node_preds = node_preds,
      tree_paths = tree_paths,
      node_ids = inner_tree_info$nodeID,
      split_vars = inner_tree_info$idx,
      split_vals = inner_tree_info$splitval
    ) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      setNames(paste0(".node", inner_tree_info$nodeID))
  } else {
    psi <- extract_psi_chr_cpp(
      x = x_mat,
      node_preds = node_preds,
      tree_paths = tree_paths,
      node_ids = inner_tree_info$nodeID,
      split_vars = inner_tree_info$idx,
      split_vals = inner_tree_info$splitval,
      unordered_factors = as.integer(
        which(colnames(x) %in% unordered_factors) - 1
      )
    ) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      setNames(paste0(".node", inner_tree_info$nodeID))
  }

  if (normalize) {
    if (is.null(psi_train)) {
      if (is.null(inbag_counts)) {
        inbag_counts <- rep(1, nrow(psi))
      }
      unique_psi_values <- purrr::map(
        psi,
        ~ c(-sqrt(sum((.x == 1) * inbag_counts) /
                    sum((.x == -1) * inbag_counts)),
            sqrt(sum((.x == -1) * inbag_counts) /
                   sum((.x == 1) * inbag_counts)))
      )
    } else {
      unique_psi_values <- purrr::map(psi_train, ~ sort(setdiff(unique(.x), 0)))
    }
    psi <- purrr::map2(
      .x = psi, .y = unique_psi_values,
      ~ dplyr::case_when(
        .x == -1 ~ .y[1],
        .x == 1 ~ .y[2],
        TRUE ~ 0
      )
    ) %>%
      as.data.frame()
  }

  return(psi)
}

#' Extract Psi (i.e., stump) matrix from tree (non-Rcpp version).
#'
#' @inheritParams extract_psi_fast
#' @export
extract_psi <- function(x, tree_info, tree_paths, node_preds, psi_train = NULL,
                        inbag_counts = NULL, normalize = FALSE,
                        unordered_factors = NULL) {
  inner_tree_info <- tree_info %>%
    dplyr::filter(!terminal)

  if (length(unordered_factors) == 0) {
    psi <- purrr::map(
      1:nrow(x),
      function(i) {
        terminal_node <- node_preds[i]
        tree_path <- tree_paths[[as.character(terminal_node)]]
        in_node <- which(inner_tree_info$nodeID %in% tree_path)
        split_dir <- c(x[i, inner_tree_info$splitvarName[in_node]] > inner_tree_info$splitval[in_node])
        psi_vec <- rep(0, nrow(inner_tree_info))
        psi_vec[in_node] <- split_dir * 2 - 1
        return(data.frame(t(psi_vec)))
      }
    )
  } else {
    psi <- purrr::map(
      1:nrow(x),
      function(i) {
        terminal_node <- node_preds[i]
        tree_path <- tree_paths[[as.character(terminal_node)]]
        in_node <- which(inner_tree_info$nodeID %in% tree_path)
        split_dir <- purrr::map2_lgl(
          inner_tree_info$splitvarName[in_node],
          inner_tree_info$splitval[in_node],
          function(split_var, split_val) {
            if (split_var %in% unordered_factors) {
              x[i, split_var] %in%
                as.numeric(stringr::str_split(split_val, ",")[[1]])
            } else {
              x[i, split_var] > as.numeric(split_val)
            }
          }
        )
        psi_vec <- rep(0, nrow(inner_tree_info))
        psi_vec[in_node] <- split_dir * 2 - 1
        return(data.frame(t(psi_vec)))
      }
    )
  }
  psi <- psi %>%
    purrr::list_rbind() %>%
    setNames(paste0(".node", inner_tree_info$nodeID))

  if (normalize) {
    if (is.null(psi_train)) {
      if (is.null(inbag_counts)) {
        inbag_counts <- rep(1, nrow(psi))
      }
      unique_psi_values <- purrr::map(
        psi,
        ~ c(-sqrt(sum((.x == 1) * inbag_counts) /
                    sum((.x == -1) * inbag_counts)),
            sqrt(sum((.x == -1) * inbag_counts) /
                   sum((.x == 1) * inbag_counts)))
      )
    } else {
      unique_psi_values <- purrr::map(psi_train, ~ sort(setdiff(unique(.x), 0)))
    }
    psi <- purrr::map2(
      .x = psi, .y = unique_psi_values,
      ~ dplyr::case_when(
        .x == -1 ~ .y[1],
        .x == 1 ~ .y[2],
        TRUE ~ 0
      )
    ) %>%
      as.data.frame()
  }

  return(psi)
}


#' Get augmented x matrix (i.e., (x_s, psi)) from x and psi.
#'
#' @inheritParams extract_psi_fast
#' @inheritParams rfplus
#' @param x A data frame or matrix of features.
#' @param psi A data frame or matrix of psi values.
#' @param normalize Whether to normalize psi values so that largest stump per
#'   feature has the same standard deviation as its raw feature.
#' @param normalize_factor A numeric vector of length equal to the number of
#'   features in x. If provided, this vector will be used to normalize the
#'   features in x. If NULL, the standard deviation of each feature in x will
#'   be used.
#'
#' @export
get_augmented_x <- function(x, psi, tree_info = NULL, include_raw = TRUE,
                            normalize = FALSE, normalize_factor = NULL) {
  if (include_raw) {
    splitvars <- tree_info %>%
      dplyr::filter(!is.na(splitvarName)) %>%
      dplyr::pull(splitvarName) %>%
      unique()
    x_splitvar <- x[, splitvars, drop = FALSE]
    x_splitvar <- dummy_code(x_splitvar)$x
    if (isFALSE(normalize)) {
      x_aug <- as.matrix(cbind(x_splitvar, psi))
    } else {
      if (any(sapply(x, class) == "factor")) {
        # TODO: fix for factor/categorical features
        stop("Cannot normalize if factor variables are included in x.")
      }
      if (is.null(normalize_factor)) {
        psi_sds <- tibble::tibble(
          nodeID = colnames(psi),
          psi_sd = apply(psi, 2, sd)
        )
        x_sds <- tibble::tibble(
          var = colnames(x_splitvar),
          x_sd = apply(x_splitvar, 2, sd)
        )
        scale_factors <- tree_info %>%
          dplyr::mutate(nodeID = paste0(".node", nodeID)) %>%
          dplyr::right_join(psi_sds, by = "nodeID") %>%
          dplyr::group_by(splitvarName) %>%
          dplyr::summarise(max_psi_sd = max(psi_sd)) %>%
          dplyr::right_join(x_sds, by = c("splitvarName" = "var")) %>%
          dplyr::mutate(scale_factor = max_psi_sd / x_sd)
        order_indices <- match(scale_factors$splitvarName, splitvars)
        scale_factors <- scale_factors %>%
          dplyr::arrange(!!order_indices) %>%
          dplyr::pull(scale_factor)
      } else {
        scale_factors <- normalize_factor
      }
      x_aug <- cbind(sweep(x_splitvar, 2, scale_factors, "*"), psi) %>%
        as.matrix()
      attr(x_aug, "x_scale_factors") <- scale_factors
    }
  } else {
    x_aug <- as.matrix(psi)
  }
  return(x_aug)
}


#' Dummy code factor variables in x and x_test.
#'
#' @param x A data frame containing the training data.
#' @param x_test A data frame containing the test data.
#'
#' @export
dummy_code <- function(x, x_test = NULL) {
  if (is.data.frame(x)) {
    if (any(sapply(x, class) == "factor")) {
      x_full <- dplyr::bind_rows(x, x_test)
      x_full <- data.frame(model.matrix(~ ., x_full)[, -1])  # remove intercept
      x <- x_full[1:nrow(x), , drop = FALSE]
      if (!is.null(x_test)) {
        x_test <- x_full[(nrow(x) + 1):nrow(x_full), , drop = FALSE]
      }
    }
  }
  return(list(x = x, x_test = x_test))
}


#' Preprocess x before computing Psi matrix.
#'
#' @param x A data frame containing the data to preprocess.
#' @param x_train A data frame containing the training data used for fitting
#'   the RF+ model.
#' @param fit An object of class \code{rfplus}.
#'
#' @export
get_x_psi <- function(x, x_train, fit) {
  x_psi <- x
  unordered_factors <- colnames(x_train)[!fit$forest$is.ordered]
  factor_cols <- sapply(x_train, class) == "factor"
  if (any(factor_cols)) {
    for (col in names(factor_cols)) {
      x_psi <- x_psi %>%
        dplyr::mutate(
          {{col}} := factor(
            as.character(.data[[col]]),
            levels = fit$forest$covariate.levels[[col]]
          ) %>%
            as.numeric()
        )
    }
  }
  return(x_psi)
}


#' Wrapper to get Psi for training and test data.
#'
#' @param object An object of class \code{rfplus}.
#' @param x A data frame containing the test data.
#' @param use_rcpp A logical indicating whether to use Rcpp code.
#'
#' @export
get_train_test_psis <- function(object, x,
                                # get rid of arguments later
                                use_rcpp = TRUE) {
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

  node_train_preds <- predict(
    rf_fit, x_train, type = "terminalNodes", num.threads = 1
  )$predictions
  node_preds <- predict(
    rf_fit, x, type = "terminalNodes", num.threads = 1
  )$predictions
  if (use_rcpp) {
    forest_paths <- get_forest_paths_fast(tree_infos)
  } else {
    forest_paths <- get_forest_paths(tree_infos)
  }

  unordered_factors <- colnames(x_train)[!rf_fit$forest$is.ordered]
  x_train_psi <- get_x_psi(x_train, x_train, rf_fit)
  psis_train <- purrr::imap(
    tree_infos,
    function(tree_info, tree_idx) {
      if (use_rcpp) {
        psi <- extract_psi_fast(
          x = x_train_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_train_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      } else {
        psi <- extract_psi(
          x = x_train_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_train_preds[, tree_idx],
          inbag_counts = rf_fit$inbag.counts[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      }
    }
  )

  x_psi <- get_x_psi(x, x_train, rf_fit)
  psis <- purrr::imap(
    tree_infos,
    function(tree_info, tree_idx) {
      if (use_rcpp) {
        psi <- extract_psi_fast(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          psi_train = psis_train[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      } else {
        psi <- extract_psi(
          x = x_psi,
          tree_info = tree_info,
          tree_paths = forest_paths[[tree_idx]],
          node_preds = node_preds[, tree_idx],
          psi_train = psis_train[[tree_idx]],
          normalize = normalize_stumps,
          unordered_factors = unordered_factors
        )
      }
    }
  )
  return(list(train = psis_train, test = psis))
}


#' Compute R-squared between two vectors while returning 0 if the estimate is
#' constant.
#'
#' @inheritParams yardstick::rsq_vec
#' @export
rsq_vec_narm <- function(truth, estimate, ...) {
  if (length(unique(estimate)) == 1) {
    return(0)
  } else {
    return(yardstick::rsq_vec(truth = truth, estimate = estimate, ...))
  }
}
