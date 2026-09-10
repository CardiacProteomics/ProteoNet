#' Construct and identify subnetworks from protein interactions
#'
#' Builds a protein-protein interaction graph from a set of interactions,
#' filters edges by confidence score, removes small connected components,
#' and partitions the remaining graph into subnetworks.
#'
#' Subnetworks are defined either by Louvain community detection
#' (\code{method = "louvain"}), which may split a component into several
#' communities, or by the connected components themselves
#' (\code{method = "connected"}), in which case each component is one cluster.
#'
#' @param interactions object of protein-protein interactions with columns
#'   \code{from}, \code{to}, and \code{score}
#' @param min_cluster_size Minimum number of nodes required to retain a connected component
#' @param score_threshold Minimum interaction score to retain an edge
#' @param mapped_proteins Data frame mapping STRING IDs to protein identifiers
#' @param method Clustering method: \code{"louvain"} (default) or \code{"connected"}
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{graph}{An igraph object of the filtered interaction network}
#'   \item{communities}{A \code{communities} object describing the partition}
#'   \item{method}{The clustering method used}
#' }
#'
#' @export

construct_network <- function(interactions, min_cluster_size, score_threshold,
                              mapped_proteins,
                              method = c("louvain", "connected")) {

  method <- match.arg(method)

  interactions <- interactions[interactions$score >= score_threshold, ]

  edges <- interactions |> dplyr::select(from = from, to = to, score = score)

  interaction_graph <- tidygraph::tbl_graph(edges = edges, directed = FALSE)

  igraph::E(interaction_graph)$linestyle <- "solid"

  comp <- igraph::components(interaction_graph)

  big_comps <- which(comp$csize >= min_cluster_size)

  interaction_graph_filtered <- igraph::induced_subgraph(
    interaction_graph,
    vids = igraph::V(interaction_graph)[comp$membership %in% big_comps]
  )

  communities <- switch(
    method,
    louvain = igraph::cluster_louvain(interaction_graph_filtered),
    connected = {
      comp_filtered <- igraph::components(interaction_graph_filtered)
      igraph::make_clusters(
        interaction_graph_filtered,
        membership = comp_filtered$membership,
        algorithm  = "connected components",
        modularity = TRUE
      )
    }
  )

  igraph::V(interaction_graph_filtered)$community <- igraph::membership(communities)

  return(list(
    graph = interaction_graph_filtered,
    communities = communities,
    method = method
  ))
}
