

place_singletons <- function( df_out, interactions, mapped_proteins ){
  
  ens_singletons <- df_out[df_out$cluster=='not_in_cluster', ]$string_name
  ens_in_cluster <- df_out[!df_out$cluster=='not_in_cluster', ]$string_name
  
  
  results <- list()
  for( ens_i in ens_singletons ){
    
    ind_1 <- interactions$from == ens_i | interactions$to == ens_i
    ind_2 <- interactions$from %in% ens_in_cluster | interactions$to %in% ens_in_cluster
    
    
    interactions_filtered <- interactions[ind_1 & ind_2, ]
    ind_max <- which.max(interactions_filtered$score)
    
    
    score <- ifelse(length(ind_max)>0, interactions_filtered$score[ind_max], NA)
    connection_to <- ifelse(length(ind_max)>0, interactions_filtered$to[ind_max], NA)
    connection_from <- ifelse(length(ind_max)>0, interactions_filtered$from[ind_max], NA)
    
    connection <- ifelse(connection_to == ens_i, connection_from, connection_to)
    con_gn <- mapped_proteins[mapped_proteins$STRING_id==connection, ]$protein
    connection_gene_out <- ifelse(length(con_gn)>0, con_gn, ens_i)
    
    
    gn <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    gene_out<- ifelse(length(gn)>0, gn, ens_i)
    
    results[[ length(results) + 1 ]] <- data.frame(ens_id = ens_i, gene_out, connection, connection_gene_out, score)  
    
  }
  
  df_singleton <- do.call(rbind, results)
  
  return(na.omit(df_singleton))
}

