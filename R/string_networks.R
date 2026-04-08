



identify_interactions <- function( my_geneset, string_database_location, species, threshold ){
  
  string_db <- STRINGdb::STRINGdb$new(version="12", species=species, score_threshold=0, input_directory=string_database_location, link_data='detailed')
  mapped_proteins <- string_db$map(data.frame(protein = my_geneset), "protein", removeUnmappedRows = TRUE)
  all_interactions_in <- string_db$get_interactions(mapped_proteins$STRING_id)
  interactions <- all_interactions_in %>% distinct(from, to, .keep_all = TRUE)
  
  interactions$score <- 1 - (1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$database/1000)
  
  interactions$alternative_score <- 1 -(1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$neighborhood/1000)*
    (1-interactions$fusion/1000)*(1-interactions$cooccurence/1000)*(1-interactions$textmining/1000)
  
    map_temp <- merge(
      mapped_proteins,
      string_db$proteins,
      by.x = "STRING_id",
      by.y = "protein_external_id"
    )
    
  mapped_proteins$protein <- map_temp$preferred_name
  mapped_proteins$STRING_id <- map_temp$STRING_id

  
  return(list(interactions=interactions, string_db = string_db, mapped_proteins = mapped_proteins))
}


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



prepare_labels <- function(cat_sub, df_ora, communities){
  
  label_out <- list()
  db_out <- c()
  label_genes <- c()
  label_genes_out <- list()
  
  
  n_communities <- length(cat_sub$communities$original_genes)
  
  
  for(i in 1:n_communities){
    label <- c()
    ### Get the label with gene names
    gene_text <- stringr::str_wrap(paste(cat_sub$communities$original_genes[[i]], collapse = " "), width = 40)
    gene_single_text <- stringr::str_wrap(paste(cat_sub$communities$singleton_genes[[i]], collapse = " "), width = 40)
    
    if(length(gene_single_text)>0){
      label_genes <- c(label_genes, paste0("\n\n", gene_text, "\n SINGLE GENES: ", gene_single_text))
    }else{
      label_genes <- c(label_genes, paste0("\n\n", gene_text))
    }
    
    label_genes_out[[i]] <- label_genes
    
    for( db_j in unique(df_ora$database) ){
      
      
      for( com_i in unique(df_ora$category_cluster)){
        
        ind <- (df_ora$category_cluster == com_i) & (df_ora$selected) & (df_ora$database == db_j) & (df_ora$community == i) 
        df_sub <- df_ora[ind, ]
        if(sum(ind)>0){
        
          if(dim(df_sub)[1]>0){
            label_i <- strsplit(df_sub$functional_category, "_")[[1]]
            label_i <- label_i[2:length(label_i)]
          }else{
            label_i <- NA
          }
          
          
          label_i <- tolower(label_i)
          label_i <- stringr::str_wrap(label_i, width = 40)
          
          
          
          db_out <- c(db_out, db_j)
          label <- c(label, paste0(paste0(label_i, collapse =" "), " (", db_j, ")", "\n"))
        }
        
      }
      
    }    
    
    #labels_k <- tolower(labels_k)
    
    label_out[[i]] <- c(paste0(label, collapse=""), "\n\n")  
  } #n_communities 
  return(label_out)
}

prepare_labels_alternative <- function(layout_fr, df_ora){
  
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
            label_i <- strsplit(df_sub$functional_category, "_")[[1]]
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



get_community_representatives <- function( cat_sub, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max, folder_databases ){
  
  n_databases <- length(databases_tested)
  
  for(j in 1:n_databases){
    
    db <- get_database(databases_tested[[j]], folder_databases)
    
    
    n_communities <- length(cat_sub$communities$all_genes)
    
    for(i in 1:n_communities){
      
      
      df_ora_all <- genes_ORA(cat_sub$communities$all_genes[[i]], universe, db, ora_min, ora_max)
      df_ora <- df_ora_all[df_ora_all$pv_bh<=0.05, ]
      
      if(dim(df_ora)[1]>1){
        overlap_matrix <- get_overlap_matrix( df_ora$functional_category, db )
        df_cluster_ora <- get_sub_clusters_ORA( overlap_matrix, df_ora, threshold_mean, threshold_min )
        df_ora$category_cluster <- df_cluster_ora$cluster
        
        if( selection == "fdr"){
          df_ora <- df_ora %>%
            group_by(category_cluster) %>%
            mutate(
              selected = pv_bh == min(pv_bh, na.rm = TRUE)
            ) %>%
            ungroup()
        }
        
        if( selection == "set_size"){
          df_ora <- df_ora %>%
            group_by(category_cluster) %>%
            mutate(
              selected = pv_bh == min(pv_bh, na.rm = TRUE)
            ) %>%
            ungroup()
        }
        
        df_ora <- as.data.frame(df_ora)
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        if(i == 1){ df_save <- df_ora }else{ df_save <- rbind(df_save, df_ora)}
      }else{
        df_ora$category_cluster <- rep(1, dim(df_ora)[1])
        df_ora$selected <- rep(TRUE, dim(df_ora)[1])
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        
        if(i == 1){ df_save <- df_ora }else{ df_save <- rbind(df_save, df_ora)}
        
      }
      
      
    }

    
    df_save$database <- rep(strsplit(names(db)[1], '_')[[1]][1], dim(df_save)[1])
    
    if(j == 1){df_ora_combined <- df_save } else{ df_ora_combined <- rbind(df_ora_combined, df_save)}
  } # j 
  
  return(df_ora_combined)
}


get_community_representatives_alternative <- function( layout_fr, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max, folder_databases ){
  
  n_databases <- length(databases_tested)
  
  for(j in 1:n_databases){
    
    db <- get_database(databases_tested[[j]], folder_databases)
    
    
    communities_un <- unique(layout_fr$community)
    
    for(i in communities_un){
      
      print(i)
      
      genes_community <- layout_fr$gene_name[layout_fr$community == i]
      
      
      df_ora_all <- genes_ORA(genes_community, universe, db, ora_min, ora_max)
      df_ora <- df_ora_all[df_ora_all$pv_bh<=0.05, ]
      
      if(dim(df_ora)[1]>1){
        overlap_matrix <- get_overlap_matrix( df_ora$functional_category, db )
        df_cluster_ora <- get_sub_clusters_ORA( overlap_matrix, df_ora, threshold_mean, threshold_min )
        df_ora$category_cluster <- df_cluster_ora$cluster
        
        if( selection == "fdr"){
          df_ora <- df_ora %>%
            group_by(category_cluster) %>%
            mutate(
              selected = pv_bh == min(pv_bh, na.rm = TRUE)
            ) %>%
            ungroup()
        }
        
        if( selection == "set_size"){
          df_ora <- df_ora %>%
            group_by(category_cluster) %>%
            mutate(
              selected = pv_bh == min(pv_bh, na.rm = TRUE)
            ) %>%
            ungroup()
        }
        
        df_ora <- as.data.frame(df_ora)
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        if(i == 1){ df_save <- df_ora }else{ df_save <- rbind(df_save, df_ora)}
      }else{
        df_ora$category_cluster <- rep(1, dim(df_ora)[1])
        df_ora$selected <- rep(TRUE, dim(df_ora)[1])
        df_ora$community <- rep(i, dim(df_ora)[1])
        
        
        if(i == 1){ df_save <- df_ora }else{ df_save <- rbind(df_save, df_ora)}
        
      }
      
      
    }
    
    
    
    df_save$database <- rep(strsplit(names(db)[1], '_')[[1]][1], dim(df_save)[1])
    
    if(j == 1){df_ora_combined <- df_save } else{ df_ora_combined <- rbind(df_ora_combined, df_save)}
  } # j 
  
  return(df_ora_combined)
}


identify_category_subclusters <- function(communities, df_singleton, string_db,interaction_graph_filtered, singleton_threshold){
  
  genes_not_in_cluster <- c()
  #genes_in_cluster_singleton <- c()
  #genes_in_cluster_OG <- c()
  #all_genes_in_cluster <- c()
  all_genes <- list()
  original_genes <- list()
  singleton_genes <- list()
  
  for (community_id in unique(communities$membership)) {
    
    #print(paste0("community", community_id))
    # Get the protein list for the current community
    community_proteins <- igraph::V(interaction_graph_filtered)$name[igraph::V(interaction_graph_filtered)$community == community_id]
    
    
    
    df <- string_db$add_proteins_description(data.frame(STRING_id = community_proteins))
    #print(community_id)
    #print(community_id)

      df$preferred_name <- df$preferred_name
    
    
    
 
    
    
    ### add back singletons to the network
    ind_cluster_single <- (df_singleton$cluster_out_max == community_id) & (df_singleton$score_out_max>singleton_threshold)
    ind_not_in_cluster_single <- (df_singleton$cluster_out_max == community_id) & (df_singleton$score_out_max<=singleton_threshold)
    
    genes_not_in_community <- df_singleton[ind_not_in_cluster_single, ]$gene_out 
    genes_single <- df_singleton[ind_cluster_single, ]$gene_out 
    genes_in_community <- df$preferred_name
    
    genes_not_in_cluster <- c(genes_not_in_cluster, genes_not_in_community )
    
    all_genes[[community_id]] <- c(df$preferred_name, genes_single)
    original_genes[[community_id]] <-  df$preferred_name
    
    if(length(genes_single>0)){
      singleton_genes[[community_id]] <-  genes_single
    }else{
      singleton_genes[[community_id]] <-  list(NA)
    }
    
  }
  
  list_out <- list(
    not_community = list(
      genes_not_in_cluster = genes_not_in_cluster
    ),
    communities = list(
      all_genes       = all_genes, 
      original_genes  = original_genes, 
      singleton_genes = singleton_genes
    )
  )
  
  return(list_out)
}


identify_category_subclusters_alternative <- function(communities, df_singleton, df_out, string_db, interaction_graph_with_singletons){
  
  singletons_all <- df_out$gene[grepl("not_in_cluster", df_out$cluster)]
  singletons_in_network <- intersect(singletons_all, layout_fr$gene_name)
  singletons_not_in_network <- setdiff(singletons_all, layout_fr$gene_name)
  
  
  genes_not_in_cluster <- c()
  #genes_in_cluster_singleton <- c()
  #genes_in_cluster_OG <- c()
  #all_genes_in_cluster <- c()
  all_genes <- list()
  original_genes <- list()
  singleton_genes <- list()
  
  for (community_id in unique(communities$membership)) {
    
    #print(paste0("community", community_id))
    # Get the protein list for the current community
    community_proteins <- igraph::V(interaction_graph_with_singletons)$name[igraph::V(interaction_graph_with_singletons)$community == community_id]
    
    
    
    df <- string_db$add_proteins_description(data.frame(STRING_id = community_proteins))
    #print(community_id)
    #print(community_id)
    
    df$preferred_name <- df$preferred_name
    
    
    genes_not_in_community <- df_singleton[ind_not_in_cluster_single, ]$gene_out 
    genes_single <- df_singleton[ind_cluster_single, ]$gene_out 
    genes_in_community <- df$preferred_name
    
    genes_not_in_cluster <- c(genes_not_in_cluster, genes_not_in_community )
    
    all_genes[[community_id]] <- df$preferred_name
    original_genes[[community_id]] <-  df$preferred_name
    
    if(length(genes_single>0)){
      singleton_genes[[community_id]] <-  genes_single
    }else{
      singleton_genes[[community_id]] <-  list(NA)
    }
    
  }
  
  list_out <- list(
    not_community = list(
      genes_not_in_cluster = genes_not_in_cluster
    ),
    communities = list(
      all_genes       = all_genes, 
      original_genes  = original_genes, 
      singleton_genes = singleton_genes
    )
  )
  
  return(list_out)
}




assign_subnetwork_membership <- function(my_geneset, interactions, communities, string_db, mapped_proteins){
  
  
  
  
  all_stringID <- c(interactions$from, interactions$to) 
  all_stringID_unique <- unique(all_stringID)
  
  genes_in_network <- c()
  network_membership <- c()
  
  
  for( i in 1:length(all_stringID_unique)){
    
    si <- communities$names[i]
    
    membership <- communities$membership[i]
    
    genes_in_network <- c(genes_in_network, mapped_proteins[mapped_proteins$STRING_id==si, ]$protein)
    network_membership <- c(network_membership, membership)
  }
  genes_not_in_network <- setdiff(my_geneset, genes_in_network)
  
  names <- c()
  cluster <- c()
  string_name <- c()
  
  ens2membership <- setNames( communities$membership, communities$names )
  
  #### make dataframe output specifying subcluster membership
  for(gene in my_geneset){
    
    ens <- mapped_proteins[mapped_proteins$protein==gene, ]$STRING_id
    
    if(length(ens)==0){
      cluster <- c(cluster, "not_in_cluster")
      string_name <- c(string_name, paste(gene, "(not_in_string)"))
      names <- c(names, gene)
    }else{
      for(ens_i in ens){
        
        if(ens_i %in% communities$names){
          
          
          membership <- ens2membership[[ens_i]]
        
          
          
          
          cluster <- c(cluster, membership)
          names <- c(names, gene)
          string_name <- c(string_name, ens_i)
          
          #print(gene)
        }else if(gene %in% genes_not_in_network){
          
          cluster <- c(cluster, "not_in_cluster")
          names <- c(names, gene)
          string_name <- c(string_name, ens_i)
          
        }else{
          #cluster <- c(cluster, "not_located")
          #print(gene)
        }
        
      }
    } 
    
  }
  
  df_out <- data.frame(cluster, gene = names, string_name )
  df_out <- df_out[order(df_out$cluster), ]
  
  return(list(df_out = df_out, mapped_proteins = mapped_proteins))
  
}


place_singletons <- function( df_out, mapped_proteins, interactions ){
  
  ens_out <- c()
  cluster_out <- c()
  score_out <- c()
  gene_out <- c()
  
  cluster_out_max <- c()
  score_out_max <- c()
  
  cluster_out_mean <- c()
  score_out_mean <- c()
  
  for(ens_i in df_out[df_out$cluster=='not_in_cluster', ]$string_name){
    
    
    ind_from <- which(interactions$from == ens_i)
    ind_to <- which(interactions$to == ens_i)
    
    from_connections <- interactions$to[ind_from]
    from_scores <- interactions$alternative_score[ind_from]
    
    to_connections <- interactions$from[ind_to]
    to_scores <- interactions$alternative_score[ind_to]
    
    
    #print(interactions[which(interactions$from == ens_i | interactions$to == ens_i), ])
    
    clusters <- unique(df_out$cluster)
    n_clusters <- length(clusters)
    
    cluster_score <- setNames(rep(0,n_clusters), clusters)
    cluster_score_max <- setNames(rep(0,n_clusters), clusters)
    cluster_score_mean <- setNames(rep(0,n_clusters), clusters)
    
    cluster_score_from <- setNames(rep(0,n_clusters), clusters) 
    
    #if(F){
    if(length(from_connections)>0){
      for( j in 1:length(from_connections)){
        
        ens_j <- from_connections[j]
        score_j <- from_scores[j]
        
        ind <- ens_j==df_out$string_name
        if(sum(ind)>0){
          
          cluster_i <- df_out$cluster[ind]
          
          cluster_score_from[[cluster_i]]<- cluster_score_from[[cluster_i]] +  score_j
          cluster_score_max[[cluster_i]]<- max(c(cluster_score_max[[cluster_i]], score_j))
          cluster_score_mean[[cluster_i]]<- cluster_score_mean[[cluster_i]] +  score_j/(length(from_connections) + length(to_connections))
          
          #print(paste(ens_j, score_j, df_out$string_name[ind], df_out$cluster[ind]))
          
        }
      }
    }
    #}
    if(length(to_connections)>0){
      cluster_score_to <- setNames(rep(0,n_clusters), clusters) 
      for( j in 1:length(to_connections)){
        ens_j <- to_connections[j]
        score_j <- to_scores[j]
        
        ind <- ens_j==df_out$string_name
        if(sum(ind)>0){
          
          cluster_i <- df_out$cluster[ind]
          n_cluster <- sum(df_out$cluster==cluster_i, na.rm=T)
          cluster_score_to[[cluster_i]]<- cluster_score_to[[cluster_i]] +  score_j
          cluster_score[[cluster_i]]<- cluster_score[[cluster_i]] +  score_j
          cluster_score_max[[cluster_i]]<- max(c(cluster_score_max[[cluster_i]], score_j))
          cluster_score_mean[[cluster_i]]<- cluster_score_mean[[cluster_i]] +  score_j/(length(from_connections) + length(to_connections))
          #print(paste(ens_i, ens_j, score_j, df_out$string_name[ind], df_out$cluster[ind]))
          
        }
      }
    }
    
    
    exclude_name <- "not_in_cluster"  
    filtered_list <- cluster_score[setdiff(names(cluster_score), exclude_name)]
    filtered_list_max <- cluster_score_max[setdiff(names(cluster_score_max), exclude_name)]  
    filtered_list_mean <- cluster_score_mean[setdiff(names(cluster_score_mean), exclude_name)]  
    
    max_name <- names(filtered_list)[which.max(unlist(filtered_list))]
    max_value <- max(unlist(filtered_list))
    
    maxmax_name <- names(filtered_list_max)[which.max(unlist(filtered_list_max))]
    maxmax_value <- max(unlist(filtered_list_max))
    
    meanmax_name <- names(filtered_list_mean)[which.max(unlist(filtered_list_mean))]
    meanmax_value <- max(unlist(filtered_list_mean))
    
    ens_out <- c(ens_out, ens_i)
    cluster_out <- c(cluster_out, max_name)
    score_out <- c(score_out, max_value)
    
    cluster_out_max <- c(cluster_out_max, maxmax_name)
    score_out_max <- c(score_out_max, maxmax_value)
    
    cluster_out_mean <- c(cluster_out_mean, meanmax_name)
    score_out_mean <- c(score_out_mean, meanmax_value)
    
    gn <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    
    gene_out<- c(gene_out, ifelse(length(gn)>0, gn, ens_i))
    
  }
  
  df_singleton <- data.frame(gene_out, cluster_out, score_out, cluster_out_max, score_out_max, cluster_out_mean, score_out_mean)
  
  return(df_singleton)
  
}

place_singletons_alternative <- function( df_out, interactions, mapped_proteins ){
  
  to_connection <- c()
  to_score <- c()
  
  from_connection <- c()
  from_score <- c()
  
  gene_out <- c()
  
  gene_to <- c()
  gene_from <- c()
  
  ens_singletons <- df_out[df_out$cluster=='not_in_cluster', ]$string_name
  ens_in_cluster <- df_out[!df_out$cluster=='not_in_cluster', ]$string_name
  
  
  
  for( ens_i in ens_singletons ){
    
    interactions_filtered <- interactions[interactions$from %in% ens_in_cluster |
                                            interactions$to %in% ens_in_cluster, ]
    
    ind_from <- which(interactions_filtered$from == ens_i)
    ind_to <- which(interactions_filtered$to == ens_i)
    
    from_connections <- interactions_filtered$to[ind_from]
    from_scores <- interactions_filtered$alternative_score[ind_from]
    
    to_connections <- interactions_filtered$from[ind_to]
    to_scores <- interactions_filtered$alternative_score[ind_to]
    
    ind_max_to <- which.max(to_scores)
    ind_max_from <- which.max(from_scores)
    
    to_score <- c(to_score, ifelse(length(ind_max_to)>0, to_scores[ind_max_to], NA))
    to_connection <- c(to_connection, ifelse(length(ind_max_to)>0, to_connections[ind_max_to], NA))
    
    from_score <- c(from_score, ifelse(length(ind_max_from)>0, from_scores[ind_max_from], NA))
    from_connection <- c(from_connection, ifelse(length(ind_max_from)>0, from_connections[ind_max_from], NA))
    
    gn <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    
    gene_out<- c(gene_out, ifelse(length(gn)>0, gn, ens_i))
    
  }
  
  
  for(i in 1:length(from_connection)){
    
    ens_to <- to_connection[i]
    ens_from <- from_connection[i]
    
    gene_to <- c(gene_to, ifelse(!ens_to=="NA", mapped_proteins[mapped_proteins$STRING_id==ens_to, ]$protein, NA))
    gene_from <- c(gene_from, ifelse(!ens_from=="NA", mapped_proteins[mapped_proteins$STRING_id==ens_from, ]$protein, NA))
    
  }
  
  df_tmp <- data.frame( ens_singletons, to_connection, from_connection, gene_out, gene_from, gene_to, from_score, to_score)
  df_singleton <- df_tmp %>%
    mutate(
      max_score = pmax(from_score, to_score, na.rm = TRUE),
      max_gene = case_when(
        !is.na(from_score) &  from_score == max_score ~ from_connection,
        !is.na(to_score)   &  to_score   == max_score ~ to_connection,
        TRUE                                          ~ NA_character_
      )
    )
  
  return(df_singleton)
}


save_network_structure <- function(cat_sub, directory, reference){
  
  gene_out <- c()
  com <- c()
  category <- c()
  
  n_communities <- length(cat_sub$communities$original_genes)
  
  original_genes <- cat_sub$communities$original_genes
  singleton_genes <- cat_sub$communities$singleton_genes
  outside_genes <- cat_sub$not_community$genes_not_in_cluster
  
  for(i in 1:n_communities){
    
    community_org_i <- original_genes[[i]]
    community_single_i <- singleton_genes[[i]]
    
    for( gene_i in community_org_i){
      
      gene_out <- c(gene_out, gene_i)
      com <- c(com, i)
      category <- c(category, "original")
      
    } # original genes
    
    for( gene_i in community_single_i){
      
      gene_out <- c(gene_out, gene_i)
      com <- c(com, i)
      category <- c(category, "singleton")
      
    } #singleton genes
  }
  
  
  
  for( gene_i in outside_genes){
    
    gene_out <- c(gene_out, gene_i)
    com <- c(com, NA)
    category <- c(category, "Not_in_network")
    
  } # singleton genes
  
  df_out <- data.frame(gene = gene_out, community = com, network_membership = category)
  
  df_out2 <- df_out[!is.na(df_out$gene), ]
  
  write.xlsx(df_out2, file = paste0(directory, "network_membership_",reference,".xlsx"))
}


add_singletons_to_graph <- function(interaction_graph_filtered, df_singleton, df_out, mapped_proteins, protein_highlight){
  
  interaction_graph_with_singletons <- interaction_graph_filtered

  #edges <- as_data_frame(interaction_graph_with_singletons, what = "edges") |>
  #  rename(
  #    source = from,
  #    target = to
  #  )
  
  #write.csv(edges, paste0(folder_results, "/edges.csv"), row.names = FALSE)
  
  
  
  df_single_filt <- df_singleton[grepl(as.character(species), df_singleton$ens_singletons) & !is.na(df_singleton$max_score), ]

  
  new_vertices <- c()
  new_community <- c()
  new_highlight <- c()
  for(ens_i in df_single_filt$ens_singletons){
    
    new_vertices <- c(new_vertices, ens_i)
    
    df_i <- df_single_filt[df_single_filt$ens_singletons == ens_i, ]
    
    new_community <- c(new_community, as.numeric(df_out[df_out$string_name == df_i$max_gene, "cluster"]))
    
    gene <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    
    new_highlight <- c(new_highlight, gene %in% protein_highlight)
    
  }

  

  
  

  interaction_graph_with_singletons <- igraph::add_vertices(
    interaction_graph_with_singletons,
    nv   = length(new_vertices),
    name = new_vertices, 
    community = new_community#, 
    #highlight = new_highlight
  )

  
  n_before <- igraph::ecount(interaction_graph_with_singletons)
  
  for( i in 1:dim(df_single_filt)[1]){
    
    interaction_graph_with_singletons <- igraph::add_edges(
      interaction_graph_with_singletons, 
      c(df_single_filt$ens_singletons[i], df_single_filt$max_gene[i]), 
      attr = list(alternative_score = df_single_filt$max_score[i])
    )
    
  }

  n_after <- igraph::ecount(interaction_graph_with_singletons)
  E(interaction_graph_with_singletons)[(n_before + 1):n_after]$linestyle <- "dashed"
  
  E(interaction_graph_with_singletons)$linestyle[is.na(E(interaction_graph_with_singletons)$linestyle)] <- "solid"
  
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
  
  
  if(length(protein_highlight)>0){
    layout_fr$highlight <- igraph::V(interaction_graph_with_singletons)$highlight
  }
  print("XXXXXXXX")
  
  return(list(layout = layout_fr, graph = interaction_graph_with_singletons))
}
