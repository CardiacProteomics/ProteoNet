#' Compute overlap matrix between enriched gene sets
#'
#' Calculates pairwise overlap between enriched pathways based on shared genes,
#' producing a matrix used for clustering and redundancy reduction of ORA results.
#'
#' Overlap can be computed using either genes observed in the input set
#' ("drawn_overlap") or full gene sets ("full_overlap"), and quantified using
#' different metrics (e.g. minimum overlap or Jaccard index).
#'
#' This corresponds to the redundancy assessment step in the ProteoNet pipeline.
#'
#' @param df_ora Data frame of enrichment results containing pathway names and
#'   gene set information (e.g. \code{geneset_drawn_in_universe})
#' @param database List of gene sets (e.g. from \code{fgsea::gmtPathways})
#' @param method Method for extracting genes for overlap calculation:
#'   \code{"drawn_overlap"} (default) or \code{"full_overlap"}
#' @param overlap Metric used to quantify overlap:
#'   \code{"min_overlap"} (default) or \code{"jaccard"}
#'
#' @return A square data frame representing the pairwise overlap matrix,
#'   with pathways as both rows and columns
#'
#' @export


get_overlap_matrix <- function( df_ora, database, method = "drawn_overlap", overlap = "min_overlap" ){
  
  pathways <- df_ora$pathway
  
  
  overlap_matrix <- matrix(0, nrow=length(pathways), ncol=length(pathways),
                           dimnames=list(names(pathways), names(pathways)))
  
  for( i in 1:length(pathways)){
    
        if(method == "drawn_overlap"){
            list_of_genes_1 <- strsplit(df_ora$geneset_drawn_in_universe[i], "/")[[1]]
        }else if(method == "full_overlap"){
            list_of_genes_1 <- database[[pathways[i]]]
            }



    
    for( j in 1:length(pathways)){
      

      if(method == "drawn_overlap"){
          list_of_genes_2 <- strsplit(df_ora$geneset_drawn_in_universe[j], "/")[[1]]
      }else if(method == "full_overlap"){
          list_of_genes_2 <- database[[pathways[j]]]
          }
      
        

        
        if(overlap == "min_overlap"){
            
            overlap_matrix[i, j] <-
                length(intersect(list_of_genes_1, list_of_genes_2))/(min(c(length(list_of_genes_1),  length(list_of_genes_2))))
           
        }else if(overlap == "jaccard"){
            overlap_matrix[i, j] <-
                length(intersect(list_of_genes_1, list_of_genes_2))/length(unique(c(list_of_genes_1, list_of_genes_2)))
        }
      
    }
  }
  
  
  M <- as.data.frame(overlap_matrix)
  
  rownames(M) <- pathways
  colnames(M) <- pathways
  
  
  return(M)
}

