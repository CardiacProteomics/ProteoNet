#' Replot and export network visualization
#'
#' Regenerates a network visualization from a precomputed graph layout,
#' applying styling options such as node highlighting, community coloring,
#' and label formatting. The function produces both a labeled and unlabeled
#' version of the network and saves them to disk.
#'
#' This function is useful for adjusting visualization parameters (e.g. label size,
#' legend size, figure dimensions) without recomputing the full ProteoNet pipeline.
#'
#' @param layout_fr A graph layout object (e.g. from \code{add_singletons_to_graph})
#'   containing node positions and attributes such as \code{community},
#'   \code{highlight}, and \code{gene_name}
#' @param folder_results Path to folder where output figures will be saved
#' @param network_height Height of the output figure (in inches)
#' @param network_width Width of the output figure (in inches)
#' @param label_size Size of node text labels
#' @param legend_size Size of legend text and title
#' @param labels Character vector of subnetwork labels for legend annotation
#'
#' @return Invisibly returns a list of ggplot objects:
#' \describe{
#'   \item{network}{Network plot without labels}
#'   \item{network_labels}{Network plot with gene labels}
#' }
#'
#' @export

  
  
replot_network <- function(layout_fr, folder_results, network_height, network_width, label_size, legend_size, labels){
  
highlight <- layout_fr$highlight
community <- layout_fr$community

if(length(highlight)>0){
  
  
  lvls <- levels(as.factor(layout_fr$community))
  keep <- labels != "\n\n" & !is.na(labels)
  
  g <- attr(layout_fr, "graph")
  
  igraph::E(g)$linestyle <- factor(
    igraph::E(g)$linestyle,
    levels = c( "solid", "dashed" )
  )
  
  attr(layout_fr, "graph") <- g
  

  
  highlight <- layout_fr$highlight
  g_figure <- ggraph::ggraph(layout_fr) +
    ggraph::geom_edge_link(edge_colour = "black", aes(edge_linetype = linestyle)) + 
    ggraph::scale_edge_linetype_manual(
      name   = "Edge type",  
      breaks = c("solid", "dashed"),
      values = c(
        solid  = "solid",
        dashed = "dashed"
      ), 
      labels = c(
        solid = "Network", 
        dashed = "Singleton"
      )
    )+
    
    # Outer ring: LMB1
    
    ggraph::geom_node_point(
      aes(fill  = as.factor(community)),
      color = "black",
      size  = 4,
      shape = 21,
      stroke = ifelse(layout_fr$highlight, 3, 0),
      #show.legend = FALSE
    ) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      legend.text  = ggplot2::element_text(size = legend_size),
      legend.title = ggplot2::element_text(size = legend_size)
    )
  
  g_figure <- g_figure + 
    ggplot2::scale_fill_discrete(
      name   = "",
      labels = labels[keep], 
      breaks = lvls[keep]
    )
  

  
  
}else{
  g_figure <- ggraph::ggraph(layout_fr)  + 
    ggraph::geom_edge_link(edge_colour = "black", aes(edge_linetype = linestyle)) +
    ggraph::geom_node_point(ggplot2::aes(color = as.factor(community)), size = 5) +
    ggplot2::theme_bw()
  
  g_figure <- g_figure + 
    ggplot2::scale_color_discrete(
      name   = "",
      labels = labels[keep], 
      breaks = lvls[keep]
    )
}

g_figure2 <- g_figure + 
  geom_node_text(
    aes(label = gene_name),
    repel = TRUE,        # avoids text overlapping
    size = label_size
  )


ggsave( plot = g_figure, filename = paste0(folder_results, "/ppi_network_replot_", reference, ".png"), width = network_width, height = network_height )
ggsave( plot = g_figure2, filename = paste0(folder_results, "/ppi_network_with_labels_replot_", reference, ".png"), width = network_width, height = network_height )


}
  




  