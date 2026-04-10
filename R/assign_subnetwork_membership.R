



assign_subnetwork_membership <- function(my_geneset, interactions, communities, string_db, mapped_proteins){

  ens2membership <- setNames( communities$membership, communities$names )
  
  results <- list()  

  for(gene in my_geneset){
    
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

