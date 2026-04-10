#' Evaluate whether a pathway cluster satisfies overlap criteria
#'
#' Determines whether a group of pathways forms a valid cluster based on
#' pairwise overlap values. A cluster is accepted if both the mean and
#' minimum overlap exceed specified thresholds.
#'
#' This function is used during hierarchical clustering to decide whether
#' a dendrogram branch should be retained as a cluster or further split.
#'
#' @param indices Integer vector of indices corresponding to pathways in the overlap matrix
#' @param M A square overlap matrix (e.g. from \code{get_overlap_matrix})
#' @param threshold_mean Minimum mean overlap required within the cluster
#' @param threshold_min Minimum pairwise overlap required within the cluster
#'
#' @return Logical value indicating whether the cluster satisfies the overlap criteria
#'
#' @export

check_cluster_correlation <- function(indices, M, threshold_mean, threshold_min) {
  
  if (length(indices) < 2) return(TRUE) # clusters of size 1 are fine
    subM <- M[indices, indices]
    avg_cor <- mean(subM[upper.tri(subM)], na.rm = TRUE)
    min_cor <- min(subM[upper.tri(subM)], na.rm = TRUE)

  return(avg_cor >= threshold_mean & min_cor>=threshold_min)

  }

