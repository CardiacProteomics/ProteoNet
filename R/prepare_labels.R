

prepare_labels<- function(layout_fr, df_ora){
  
  label_out <- list()
  db_out <- c()
  label_genes <- c()
  label_genes_out <- list()
  
  
  communities_un <- unique(layout_fr$community)
  
  
  for(i in communities_un){
    
    
    genes_community <- layout_fr$gene_name[layout_fr$community == i]
    
    label <- c()
    ### Get the label with gene names
    gene_text <- stringr::str_wrap(paste(genes_community, collapse = " "), width = 40)
    #gene_single_text <- stringr::str_wrap(paste(cat_sub$communities$singleton_genes[[i]], collapse = " "), width = 40)
    

      label_genes <- c(label_genes, paste0("\n\n", gene_text))

    
    label_genes_out[[i]] <- label_genes
    
    for( db_j in unique(df_ora$database) ){
      
      
      
      
      for( com_i in unique(df_ora$category_cluster)){
        
        ind <- (df_ora$category_cluster == com_i) & (df_ora$selected) & (df_ora$database == db_j) & (df_ora$community == i) 
        df_sub <- df_ora[ind, ]
        if(sum(ind)>0){
          
          if(dim(df_sub)[1]>0){
            label_i <- strsplit(df_sub$pathway, "_")[[1]]
            label_i <- label_i[2:length(label_i)]
          }else{
            label_i <- NA
          }
          
          
          label_i <- tolower(label_i)
          
          label_i <- paste0(paste0(label_i, collapse =" "), " (", db_j, ")")
          label_i <- stringr::str_wrap(label_i, width = 40, whitespace_only = TRUE)
          label_i <- paste0(label_i, "\n")

          
          db_out <- c(db_out, db_j)
          label <- c(label, label_i)
         
   

        }
        
      }
      
    }  
    label_col <- c(paste0(label, collapse=""))
    #label_col <- stringr::str_wrap(label_col, width = 40)
    label_col <- paste0(label_col, "\n")
    print(label_col )
    #labels_k <- tolower(labels_k)
    
    label_out[[i]] <- c(label_col)  
  } #n_communities 
  return(label_out)
}


