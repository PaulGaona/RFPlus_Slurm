% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFPlus-utils.R
\name{get_forest_paths_fast}
\alias{get_forest_paths_fast}
\title{Extract all root-to-leaf paths in a forest.}
\usage{
get_forest_paths_fast(tree_infos)
}
\arguments{
\item{tree_infos}{List of size num.trees with each entry being the output of
\code{\link[ranger:treeInfo]{ranger::treeInfo()}}.}
}
\value{
A list of size num.trees with each entry being a list of decision
paths in each tree.
}
\description{
Extract all root-to-leaf paths in a forest.
}
