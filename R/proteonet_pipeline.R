#' Run the ProteoNet analysis pipeline
#'
#' Executes the full ProteoNet workflow for a set of input proteins,
#' including network construction, subnetwork detection, enrichment analysis,
#' redundancy reduction, and visualization.
#'
#' The pipeline integrates STRING-based interaction networks with
#' over-representation analysis (ORA) to provide interpretable biological
#' context for significant protein hits.
#'
#' @param reference Character string used to label output files
#' @param genes_drawn Character vector of proteins/genes of interest
#' @param species Numeric species identifier (e.g. 9606 for Homo sapiens)
#' @param min_cluster_size Minimum size of subnetworks to retain
#' @param score_threshold Minimum interaction score for network edges
#' @param selection Method for selecting representative enrichment terms
#' @param databases_tested List of gene set databases to use
#' @param ora_min Minimum gene set size for ORA
#' @param ora_max Maximum gene set size for ORA
#' @param folder_databases Path to gene set database files
#' @param folder_results Path to save enrichment results
#' @param folder_figures Path to save generated figures
#'
#' @return A list containing:
#' \describe{
#'   \item{plots}{Network plots (with and without labels)}
#'   \item{enrichment}{ORA results (combined across subnetworks)}
#'   \item{network}{Graph object with subnetworks and singleton integration}
#' }
#'
#' @export



proteonet_pipeline <- function( reference,
                 genes_drawn,
                 species,
                 min_cluster_size,
                 score_threshold,
                  singleton_threshold, 
                 selection,
                 databases_tested,
                 ora_min,
                 ora_max,
                 folder_string,
                 folder_genesets,
                 folder_results,
                 folder_figures,
                 universe,
                 threshold_mean,
                 threshold_min
                 ){


  out_ii <- identify_interactions( genes_drawn,
                                   folder_string,
                                   species )

  out_cn <- construct_network( out_ii$interactions,
                               min_cluster_size,
                               score_threshold,
                               out_ii$mapped_proteins )

  df_asm <- assign_subnetwork_membership(genes_drawn,
                                         out_ii$interactions,
                                         out_cn$communities,
                                         out_ii$mapped_proteins)

  df_ps <- place_singletons(  df_asm,
                              out_ii$interactions,
                              out_ii$mapped_proteins, 
                              singleton_treshold  )

  out_astg <- add_singletons_to_graph(out_cn$graph,
                                      df_ps,
                                      df_asm,
                                      out_ii$mapped_proteins)


  df_gcr <- get_community_representatives( out_astg$layout,
                                           universe,
                                           selection,
                                           databases_tested,
                                           threshold_mean,
                                           threshold_min,
                                           ora_min,
                                           ora_max,
                                           folder_genesets)


  write.csv(df_gcr, file = paste0(folder_results, "/overrepresentation_analysis_alt_", reference, ".csv"))
  labels <- prepare_labels(out_astg$layout, df_gcr)

  plots <- produce_network_plot( out_astg$layout, labels)

  ggplot2::ggsave(plots$network_figure, file = paste0(folder_figures, "/network_figure_", reference, ".png"))
  ggplot2::ggsave(plots$network_figure_labels, file = paste0(folder_figures, "/network_figure_labels_", reference, ".png"))

  return(list(ORA_analysis = df_gcr, network = out_astg, labels = labels))
}



