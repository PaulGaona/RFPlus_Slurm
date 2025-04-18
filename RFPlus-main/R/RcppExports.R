# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

extract_psi_cpp <- function(x, node_preds, tree_paths, node_ids, split_vars, split_vals) {
    .Call('_RFPlus_extract_psi_cpp', PACKAGE = 'RFPlus', x, node_preds, tree_paths, node_ids, split_vars, split_vals)
}

extract_psi_chr_cpp <- function(x, node_preds, tree_paths, node_ids, split_vars, split_vals, unordered_factors) {
    .Call('_RFPlus_extract_psi_chr_cpp', PACKAGE = 'RFPlus', x, node_preds, tree_paths, node_ids, split_vars, split_vals, unordered_factors)
}

get_tree_paths_cpp <- function(terminal_node_ids, left_child_ids, right_child_ids, node_ids) {
    .Call('_RFPlus_get_tree_paths_cpp', PACKAGE = 'RFPlus', terminal_node_ids, left_child_ids, right_child_ids, node_ids)
}

