
check_cluster_correlation <- function(indices, M, threshold_mean, threshold_min) {
  
  if (length(indices) < 2) return(TRUE) # clusters of size 1 are fine
    subM <- M[indices, indices]
    avg_cor <- mean(subM[upper.tri(subM)], na.rm = TRUE)
    min_cor <- min(subM[upper.tri(subM)], na.rm = TRUE)

  return(avg_cor >= threshold_mean & min_cor>=threshold_min)

  }

