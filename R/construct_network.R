

construct_network <- function( interactions, min_cluster_size, string_db, protein_highlight, score_threshold, mapped_proteins){
  
  interactions <- interactions[interactions$score>=score_threshold, ]
  
  edges <- interactions %>% dplyr::select(from = from, to = to, score = score)
  
  interaction_graph <- tidygraph::tbl_graph(edges = edges, directed = FALSE)
  
  comp <- igraph::components(interaction_graph)
  
  big_comps <- which(comp$csize >= min_cluster_size)
  
  interaction_graph_filtered <- igraph::induced_subgraph(
    interaction_graph,
    vids = igraph::V(interaction_graph)[comp$membership %in% big_comps]
  )
  

  
  
  communities <- igraph::cluster_louvain(interaction_graph_filtered)
  igraph::V(interaction_graph_filtered)$community <- communities$membership
  
  
  
  if(length(protein_highlight)>0){
    highlight <- c()
    
    for(ens_i in igraph::V(interaction_graph_filtered)$name){
      
      prot_i <- mapped_proteins[mapped_proteins$STRING_id == ens_i, ]$protein
      highlight <- c(highlight, prot_i %in% protein_highlight )
   
    }
    igraph::V(interaction_graph_filtered)$highlight <- highlight
    }
  
  if(length(protein_highlight)>0){
    return(list(
      graph = interaction_graph_filtered,
      communities = communities, 
      highlights = highlight
    ))
  }else{
    return(list(
      graph = interaction_graph_filtered,
      communities = communities
    ))
  }
}

