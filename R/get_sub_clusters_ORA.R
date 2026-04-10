#' Cluster enriched pathways based on gene set overlap
#'
#' Groups enriched pathways into clusters based on their pairwise overlap,
#' using hierarchical clustering and adaptive splitting of the dendrogram.
#'
#' Clusters are defined such that pathways within a cluster satisfy specified
#' overlap criteria, enabling reduction of redundant enrichment terms.
#'
#' This corresponds to the redundancy reduction step in the ProteoNet pipeline.
#'
#' @param M A square overlap matrix (e.g. from \code{get_overlap_matrix}),
#'   with pathways as rows and columns
#' @param df_fa_in Data frame of enrichment results containing pathway-level
#'   statistics (e.g. p-values, gene counts)
#' @param threshold_mean Minimum mean overlap required within a cluster
#' @param threshold_min Minimum pairwise overlap required within a cluster
#'
#' @return A data frame with pathway clustering results, including:
#' \describe{
#'   \item{pw}{Pathway name}
#'   \item{cluster}{Assigned cluster identifier}
#'   \item{fdr}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{set_size}{Size of the gene set}
#'   \item{set_size_in_universe}{Size of the gene set within the universe}
#'   \item{genes_in_set}{Number of overlapping genes}
#' }
#'
#' @export

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

