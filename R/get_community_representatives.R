

get_community_representatives <- function( layout_fr, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max, folder_databases ){
  
  n_databases <- length(databases_tested)
  df_ora_combined <- NULL
  
  for(j in 1:n_databases){
    
    db <- get_database(databases_tested[[j]], folder_databases)
    
    df_save <- NULL
    communities_un <- unique(layout_fr$community)
    
    for(i in communities_un){
      
      print(i)
      
      genes_community <- layout_fr$gene_name[layout_fr$community == i]
      
      
      df_ora_all <- genes_ORA(genes_community, universe, db, ora_min, ora_max)
      df_ora <- df_ora_all[df_ora_all$pv_bh<=0.05, ]
      
      if(dim(df_ora)[1]>1){
        overlap_matrix <- get_overlap_matrix_2( df_ora, db )
        df_cluster_ora <- get_sub_clusters_ORA( overlap_matrix, df_ora, threshold_mean, threshold_min )
        df_ora$category_cluster <- df_cluster_ora$cluster
        
        if( selection == "fdr"){
          df_ora <- df_ora %>%
            dplyr::group_by(category_cluster) %>%
            dplyr::arrange(pv_bh) %>%
            dplyr::mutate(selected = dplyr::row_number() == 1) %>%
            dplyr::ungroup()
        }
        
        #if( selection == "set_size"){
        #  df_ora <- df_ora %>%
        #    dplyr::group_by(category_cluster) %>%
        #    dplyr::mutate(
        #      selected = pv_bh == min(pv_bh, na.rm = TRUE)
        #    ) %>%
        #    dplyr::ungroup()
        #}
        
        df_ora <- as.data.frame(df_ora)
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        df_save <- rbind(df_save, df_ora)
      }else{
        df_ora$category_cluster <- rep(1, dim(df_ora)[1])
        df_ora$selected <- rep(TRUE, dim(df_ora)[1])
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        
        df_save <- rbind(df_save, df_ora)
        
      }
    }

    df_save$database <- rep(strsplit(names(db)[1], '_')[[1]][1], dim(df_save)[1])
    
    df_ora_combined <- rbind(df_ora_combined, df_save)
  } # j 
  
  return(df_ora_combined)
}

