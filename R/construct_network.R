#' Construct and identify subnetworks from protein interactions
#'
#' Builds a protein-protein interaction graph from a set of interactions,
#' filters edges by confidence score, removes small connected components,
#' and detects communities using the Louvain algorithm.
#'
#' This corresponds to the subnetwork detection step in the ProteoNet pipeline.
#'
#' @param interactions object of protein-protein interactions with columns
#'   \code{from}, \code{to}, and \code{score}
#' @param min_cluster_size Minimum number of nodes required to retain a connected component
#' @param score_threshold Minimum interaction score to retain an edge
#' @param mapped_proteins Data frame mapping STRING IDs to protein identifiers
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{graph}{An igraph object of the filtered interaction network}
#'   \item{communities}{A community detection object from \code{igraph::cluster_louvain}}
#' }
#'
#' @export

construct_network <- function( interactions, min_cluster_size, score_threshold, mapped_proteins){
  
  interactions <- interactions[interactions$score>=score_threshold, ]
  
  edges <- interactions |> dplyr::select(from = from, to = to, score = score)
  
  interaction_graph <- tidygraph::tbl_graph(edges = edges, directed = FALSE)
  
  comp <- igraph::components(interaction_graph)
  
  big_comps <- which(comp$csize >= min_cluster_size)
  
  interaction_graph_filtered <- igraph::induced_subgraph(
    interaction_graph,
    vids = igraph::V(interaction_graph)[comp$membership %in% big_comps]
  )
  
  communities <- igraph::cluster_louvain(interaction_graph_filtered)
  igraph::V(interaction_graph_filtered)$community <- communities$membership
  
    return(list(
      graph = interaction_graph_filtered,
      communities = communities
    ))
    
}

