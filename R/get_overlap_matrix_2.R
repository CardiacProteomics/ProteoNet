

get_overlap_matrix_2 <- function( df_ora, database ){
  
  n <- dim(df_ora)[1]
  pathways <- df_ora$pathway
  
  overlap_matrix <- matrix(0, nrow=n, ncol=n, 
                           dimnames=list(pathways, pathways))
  
  for( i in 1:n){
    
    
    list_of_genes_1 <- strsplit(df_ora$geneset_drawn_in_universe[i], "/")[[1]]
    
    go_term_1 <- pathways[i]
    
    
    for( j in 1:n){
      
      
      list_of_genes_2 <- strsplit(df_ora$geneset_drawn_in_universe[j], "/")[[1]]
      
      
      go_term_2 <- pathways[j]
      
      overlap_matrix[i, j] <- length(intersect(list_of_genes_1, list_of_genes_2))/(min(c(length(list_of_genes_1),  length(list_of_genes_2))))
      
    }
  }
  
  
  M <- as.data.frame(overlap_matrix)
  
  rownames(M) <- pathways
  colnames(M) <- pathways
  
  
  return(M)
}

