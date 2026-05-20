#' Create network visualization of ProteoNet results
#'
#' Generates network plots of the protein interaction graph, with nodes colored
#' by subnetwork (community) and edges styled to distinguish original interactions
#' from singleton connections.
#'
#' Optionally adds gene labels to nodes for improved interpretability.
#'
#' @param layout_fr A graph layout object (e.g. from \code{add_singletons_to_graph})
#'   containing node positions and annotations, including \code{community},
#'   \code{gene_name}, and edge attributes
#' @param labels Character vector of labels for subnetworks (e.g. from
#'   \code{prepare_labels}), used in the legend
#'
#' @return A list with two ggplot objects:
#' \describe{
#'   \item{network_figure}{Network plot without node labels}
#'   \item{network_figure_labels}{Network plot with gene labels added}
#' }
#'
#' @export


produce_network_plot <- function( layout_fr, labels ){

  edge_df <- ggraph::get_edges()(layout_fr)

  has_linestyle <- "linestyle" %in% colnames(edge_df)

  if(has_linestyle){

    edge_geom <- ggraph::geom_edge_link(
      edge_colour = "black",
      ggplot2::aes(edge_linetype = linestyle)
    )

  } else {

    edge_geom <- ggraph::geom_edge_link(
      edge_colour = "black",
      linetype = "solid"
    )

  }
  
  g_figure <- ggraph::ggraph(layout_fr)  + 
    edge_geom +
    ggraph::geom_node_point(ggplot2::aes(color = as.factor(community)), size = 5) +
    ggplot2::theme_bw()+ 
    ggraph::scale_edge_linetype_manual(
      name = "Edge type",
      values = c(
        "solid" = "solid",
        "dashed" = "dashed"
      ), 
      labels = c(
        "solid" = "Original network",
        "dashed" = "Singleton"
        
      ))


  g_figure <- g_figure + 
    ggplot2::scale_color_discrete(
      name   = "",
      labels = labels
    )
  
  g_figure2 <- g_figure + 
    ggraph::geom_node_text(
      ggplot2::aes(label = gene_name),
      repel = TRUE,        # avoids text overlapping
      size = 5
    ) + 
    ggplot2::theme_bw()

  return(list(network_figure = g_figure, network_figure_labels = g_figure2))

}
