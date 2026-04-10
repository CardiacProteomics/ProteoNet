
library(dplyr)
library(ggraph)



automated_functional_analysis <- function( my_geneset, 
                               universe, 
                               threshold_mean, 
                               threshold_min, 
                               databases_tested, 
                               min_cluster, 
                               species, 
                               ora_min, 
                               ora_max, 
                               folder_databases,
                               folder_results,
                               reference, 
                               string_score_threshold,
                               singleton_threshold,
                               network_height, 
                               network_width,
                               protein_highlight = c()  ){
  

 
  
  string_database_location <- paste0(folder_databases, "/String")
  

  out_ii <- identify_interactions( my_geneset, string_database_location, species )  
  
  mapped_proteins <- out_ii$mapped_proteins
  interactions <- out_ii$interactions
  string_db <- out_ii$string_db
  
  #mapped_proteins <- string_db$map(data.frame(protein = my_geneset), "protein", removeUnmappedRows = TRUE)


  out_cn <- construct_network( interactions, min_cluster, string_db, protein_highlight, string_score_threshold, mapped_proteins )

  communities <- out_cn$communities
  interaction_graph_filtered <- out_cn$graph
  
  # Find singletons 
  all_stringID <- c(interactions$from, interactions$to) 
  all_stringID_unique <- unique(all_stringID)

  genes_in_network <- c()
  network_membership <- c()
  
  df_out <- assign_subnetwork_membership(my_geneset, interactions, communities, string_db, mapped_proteins)
  
  if((sum(df_out$cluster=="not_in_cluster")<(length(df_out$cluster))-1) ){
  
    df_singleton <- place_singletons( df_out, interactions, mapped_proteins  )

  
  selection <- "fdr"

  #cat_sub <- identify_category_subclusters(communities, df_singleton, string_db, interaction_graph_filtered,  singleton_threshold)
  
  res <- add_singletons_to_graph(interaction_graph_filtered, df_singleton, df_out,  mapped_proteins, protein_highlight )
  layout_fr <- res$layout
  interaction_graph_with_singletons <- res$graph

  df_ora <- get_community_representatives( layout_fr, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max, folder_databases )
  labels <- prepare_labels_alternative(layout_fr, df_ora)
  
  write.csv(df_ora, file = paste0(folder_results, "/overrepresentation_analysis_alt_", reference, ".csv"))
  
  #save_network_structure(cat_sub, folder_results, reference)


  if(length(protein_highlight)>0){
    g_figure <- ggraph::ggraph(layout_fr) +
      ggraph::geom_edge_link(edge_colour = "black", aes(edge_linetype = linestyle)) +
      
      # Outer ring: LMB1
      
      ggraph::geom_node_point(
        aes(fill  = as.factor(community)),
        color = "black",
        size  = 4,
        shape = 21,
        stroke = ifelse(layout_fr$highlight, 3, 0),
        #show.legend = FALSE
      ) + 
      ggplot2::theme_bw()
    
    
    #lvls <- levels(as.factor(layout_fr$community))
    #keep <- labels != "\n\n" & !is.na(labels)
    
    g_figure <- g_figure + 
      ggplot2::scale_fill_discrete(
        name   = "",
        labels = labels
      )
    
    
  }else{
  g_figure <- ggraph::ggraph(layout_fr)  + 
    ggraph::geom_edge_link(edge_colour = "black", aes(edge_linetype = linestyle)) +
    ggraph::geom_node_point(ggplot2::aes(color = as.factor(community)), size = 5) +
    ggplot2::theme_bw()
  
  
  
  g_figure <- g_figure + 
    ggplot2::scale_color_discrete(
      name   = "",
      labels = labels
    )
  }
  
  
  g_figure2 <- g_figure + 
    geom_node_text(
      aes(label = gene_name),
      repel = TRUE,        # avoids text overlapping
      size = 5
    ) + 
    theme_bw() #+
    if(F){
    theme(
      legend.text = ggtext::element_markdown(),  # <- this enables HTML/Markdown rendering
      axis.text.x = ggtext::element_markdown(), 
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(), 
      legend.position = "none"
    ) 
    }
  
  ggsave( plot = g_figure, filename = paste0(folder_results, "/ppi_network_", reference, ".png"), width = network_width, height = network_height )
  ggsave( plot = g_figure2, filename = paste0(folder_results, "/ppi_network_with_labels_", reference, ".png"), width = network_width, height = network_height )
  
  return(list(figure = g_figure, results = df_ora, network_layout = layout_fr, labels = labels, interactions_specifications = interaction_graph_with_singletons, n_success = TRUE))
  
  }else{
    
  return(list(n_success = FALSE))
    
  }
  
}

  

  
  
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
      name   = "Edge type",   # ← new legend title
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
  




  