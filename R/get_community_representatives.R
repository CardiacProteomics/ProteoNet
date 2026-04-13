#' Identify representative enrichment terms for subnetworks
#'
#' Performs over-representation analysis (ORA) for each subnetwork (community)
#' and identifies representative biological terms by clustering redundant terms
#' based on gene set overlap.
#'
#' For each cluster of related terms, a representative term is selected
#' (e.g. lowest adjusted p-value), enabling concise interpretation of
#' subnetwork-specific biological functions.
#'
#' This corresponds to the enrichment and redundancy reduction step
#' in the ProteoNet pipeline.
#'
#' @param layout_fr A layout object containing node annotations, including
#'   \code{community} and \code{gene_name}
#' @param universe Character vector of background genes/proteins
#' @param selection Method for selecting representative terms within clusters
#'   (e.g. \code{"fdr"})
#' @param databases_tested List of database identifiers to test
#' @param threshold_mean Mean overlap threshold for clustering enrichment terms
#' @param threshold_min Minimum overlap threshold for clustering
#' @param ora_min Minimum gene set size for ORA
#' @param ora_max Maximum gene set size for ORA
#' @param folder_databases Path to folder containing annotation databases
#'
#' @return A data frame combining enrichment results across all communities
#'   and databases, including:
#' \describe{
#'   \item{community}{Subnetwork identifier}
#'   \item{database}{Annotation database used}
#'   \item{category_cluster}{Cluster of related enrichment terms}
#'   \item{selected}{Logical indicator for representative terms}
#'   \item{...}{Additional ORA output columns (e.g. p-values, gene sets)}
#' }
#'
#' @export

get_community_representatives <- function( layout_fr, universe, selection, databases_tested, threshold_mean, threshold_min, ora_min, ora_max, folder_genesets,
    redundancy_method = "drawn_overlap",
    overlap_method = "min_overlap"){

  n_databases <- length(databases_tested)
  df_ora_combined <- NULL

  for(j in 1:n_databases){

    db <- get_database(databases_tested[[j]], folder_genesets)

    df_save <- NULL
    communities_un <- unique(layout_fr$community)

    for(i in communities_un){

      print(i)

      genes_community <- layout_fr$gene_name[layout_fr$community == i]


      df_ora_all <- genes_ORA(genes_community, universe, db, ora_min, ora_max)
      df_ora <- df_ora_all[df_ora_all$pv_bh<=0.05, ]

      if(dim(df_ora)[1]>1){
        overlap_matrix <- get_overlap_matrix( df_ora, db, redundancy_method, overlap_method )
        df_cluster_ora <- get_sub_clusters_ORA( overlap_matrix, df_ora, threshold_mean, threshold_min )
        df_ora$category_cluster <- df_cluster_ora$cluster

        if( selection == "fdr"){
          df_ora <- df_ora |>
            dplyr::group_by(category_cluster) |>
            dplyr::arrange(pv_bh) |>
            dplyr::mutate(selected = dplyr::row_number() == 1) |>
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

