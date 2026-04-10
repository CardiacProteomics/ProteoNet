
add_singletons_to_graph <- function(interaction_graph_filtered, df_singleton, df_out, mapped_proteins, protein_highlight){
  
  df_out <- df_out[!is.na(df_out$gene), ]
  
  new_vertices <- c()
  new_community <- c()
  new_highlight <- c()
  for(ens_i in df_singleton$ens_id){
    
    new_vertices <- c(new_vertices, ens_i)
    
    df_i <- df_singleton[df_singleton$ens_id == ens_i, ]
    
    new_community <- c(new_community, as.numeric(df_out[df_out$string_name == df_i$connection, "cluster"]))
    
    gene <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    
    
  }

  interaction_graph_with_singletons <- igraph::add_vertices(
    interaction_graph_filtered,
    nv   = length(new_vertices),
    name = new_vertices, 
    community = new_community
  )
  

  for( i in 1:dim(df_singleton)[1]){
    
    print(paste(df_singleton$ens_id[i], df_singleton$connection[i]))
    
    interaction_graph_with_singletons <- igraph::add_edges(
      interaction_graph_with_singletons, 
      c(df_singleton$ens_id[i], df_singleton$connection[i]), 
      attr = list(alternative_score = df_singleton$score[i], 
                  linestyle  = "dashed" )
    )
    
  }

  
  E(interaction_graph_with_singletons)$linestyle[
    is.na(E(interaction_graph_with_singletons)$linestyle)
  ] <- "solid"
  
  layout_fr <- ggraph::create_layout(
    interaction_graph_with_singletons,
    layout = "fr"
  )
  
  
  gene_name <- c()
  for(str_i in layout_fr$name){
    print(paste(str_i, mapped_proteins$protein[mapped_proteins$STRING_id == str_i]))
    gene_name <- c(gene_name, mapped_proteins$protein[mapped_proteins$STRING_id == str_i])
  }
  
  layout_fr$gene_name <- gene_name
  

  
  return(list(layout = layout_fr, graph = interaction_graph_with_singletons))
}
