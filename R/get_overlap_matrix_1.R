


get_overlap_matrix_1 <- function( pws, database ){
  

  overlap_matrix <- matrix(0, nrow=length(pws), ncol=length(pws), 
                           dimnames=list(names(pws), names(pws)))
  
  for( i in 1:length(pws)){
    
 
      list_of_genes_1 <- database[[pws[i]]]

    go_term_1 <- pws[i]

    
    for( j in 1:length(pws)){
      

        list_of_genes_2 <- database[[pws[j]]]
      
        
        go_term_2 <- pws[j]
      overlap_matrix[i, j] <- length(intersect(list_of_genes_1, list_of_genes_2))/(min(c(length(list_of_genes_1),  length(list_of_genes_2))))
      
    }
  }
  
  
  M <- as.data.frame(overlap_matrix)
  
  rownames(M) <- pws
  colnames(M) <- pws
  
  
  return(M)
}

