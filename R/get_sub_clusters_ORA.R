

get_sub_clusters_ORA <- function( M, df_fa_in, threshold_mean, threshold_min ){
  
  
  genes_in_set <- df_fa_in$gene_count
  set_size <- df_fa_in$set_size
  set_size_in_universe <- df_fa_in$set_size_in_universe
  fdr <- df_fa_in$pv_bh
  pws <- df_fa_in$pathway
  hc <- hclust(d = as.dist(1-M), method = "ward.D2")
  dend <- as.dendrogram(hc)
  leaves <- labels(dend)
  idx <- match(leaves, rownames(M))
  
  count <- 1
  clusters <- list()
  
  dend <- as.dendrogram(hc)
  not_yet_clusters <- list(dend)
  
  
  
  while (length(not_yet_clusters) > 0) {
    
    # Take the first element from the queue
    
    nyc <- not_yet_clusters[[1]]
    not_yet_clusters <- not_yet_clusters[-1]
    
    # Get indices of leaves under this branch
    leaves <- labels(nyc)
    idx <- match(leaves, rownames(M))
    
    
    # Check correlation or stop if leaf
    if (check_cluster_correlation(idx, M, threshold_mean , threshold_min ) || is.leaf(nyc)) {
      clusters[[count]] <- idx   # store indices (or nyc if you want sub-dendrograms)
      count <- count + 1
    } else {
      # Add children back to queue for further checking
      not_yet_clusters <- c(not_yet_clusters, list(nyc[[1]], nyc[[2]]))
    }
  }
  
  cluster_labels <- rep(NA, nrow(M))
  
  for (i in seq_along(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }
  
  
  
  ord <- hc$order  # keep row order consistent with heatmap
  
  
  
  term_clusters <- data.frame(
    pw        = pws[ord],
    cluster   = cluster_labels[ord],
    fdr       = fdr[ord],
    set_size  = set_size[ord], 
    set_size_in_universe  = set_size_in_universe[ord], 
    genes_in_set = genes_in_set[ord]
  )

  
  return(term_clusters)
}

