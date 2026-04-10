#' Assign genes to subnetwork communities
#'
#' Maps a set of genes or proteins to STRING identifiers and assigns each
#' to a subnetwork (community) based on precomputed community detection results.
#'
#' Genes that are not present in the STRING network or not assigned to a
#' community are labeled accordingly.
#'
#' This corresponds to the mapping step between input genes and subnetworks
#' in the ProteoNet pipeline.
#'
#' @param proteins Character vector of gene or protein identifiers
#' @param interactions Data frame of protein-protein interactions (unused but kept for pipeline consistency)
#' @param communities Community detection result (e.g. from \code{igraph::cluster_louvain})
#' @param mapped_proteins Data frame mapping gene/protein names to STRING IDs
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{cluster}{Assigned community (or "not_in_cluster")}
#'   \item{gene}{Input gene or protein identifier}
#'   \item{string_name}{Corresponding STRING identifier (or annotated missing value)}
#' }
#'
#' @export



assign_subnetwork_membership <- function(proteins, interactions, communities, mapped_proteins){

  ens2membership <- setNames( communities$membership, communities$names )
  
  results <- list()  

  for(gene in proteins){
    
    ens_ids <- mapped_proteins$STRING_id[mapped_proteins$protein==gene]
    
    if(length(ens_ids)==0){
      results[[length(results) + 1]] <- data.frame(
        cluster = "not_in_cluster", 
        gene = gene, 
        string_name = paste(gene, "(not_in_string)")
      )} else{
        
        for(ens_i in ens_ids){
          if(ens_i %in% names(ens2membership)){
            cluster <- ens2membership[[ens_i]]
          } else {
            cluster <- "not_in_cluster"
          }
          
        }
        results[[length(results) + 1]] <- data.frame(
          cluster = cluster, 
          gene = gene, 
          string_name = ens_i
        )
        
      }
  }
    
  df_out <- do.call(rbind, results)
  df_out <- df_out[order(df_out$cluster), ]
  
  return(df_out)
  
}

