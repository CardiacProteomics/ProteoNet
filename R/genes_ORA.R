
genes_ORA <- function(genes_drawn, universe, database, min_size, max_size){
  
  genes_drawn <- c(na.omit(genes_drawn))
  universe <- c(na.omit(universe))
  
  genes_drawn_in_universe <- intersect(genes_drawn, universe)
  
  N <- length( universe )
  n <- length( genes_drawn_in_universe )
  
  results <- list()
  
  i <- 1
  
  for(pw_i in names(database)){

    genes_in_pathway <- unname(database[pw_i])[[1]]
    genes_in_pathway_in_universe <- intersect(genes_in_pathway, universe)
    

    K <- length(genes_in_pathway_in_universe)
    k <- length(intersect(genes_in_pathway_in_universe, genes_drawn))
    
    if((K>=min_size) & (K<=max_size)){
      
      pvalue <- Rmpfr::asNumeric(exp(Rmpfr::mpfr(phyper(k-1, K, N-K, n, lower.tail = FALSE, log.p = TRUE), 2000)))
      
      results[[i]] <- data.frame(
        pathway = pw_i, 
        pvalue, 
        gene_ratio = as.double(k)/as.double(K), 
        geneset_drawn = paste(intersect(genes_in_pathway, genes_drawn), collapse="/"), 
        geneset_drawn_in_universe = paste(intersect(genes_in_pathway_in_universe, genes_drawn), collapse="/"), 
        gene_count = as.double(k), 
        set_size = length(genes_in_pathway), 
        set_size_in_universe = K
      )
      i <- i+1
    }  
  }

  df_out <- do.call(rbind, results)

  df_out$pv_bh <- p.adjust(df_out$pvalue, method="BH")
  
  return(df_out[order(log(df_out$pvalue)), ])
}





