#' Assign singleton proteins to nearest subnetworks
#'
#' For proteins not assigned to any subnetwork ("singletons"), identifies their
#' strongest interaction with proteins that are part of a cluster and assigns
#' them to a nearby subnetwork context.
#'
#' This improves coverage of the network by linking isolated proteins to their
#' most relevant neighbors based on interaction strength.
#'
#' @param df_out Data frame containing subnetwork membership assignments,
#'   including a \code{cluster} column and STRING identifiers
#' @param interactions Data frame of protein-protein interactions with columns
#'   \code{from}, \code{to}, and \code{score}
#' @param mapped_proteins Data frame mapping STRING IDs to gene or protein names
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{ens_id}{STRING identifier of the singleton protein}
#'   \item{gene_out}{Gene/protein name of the singleton}
#'   \item{connection}{STRING ID of the strongest connected protein}
#'   \item{connection_gene_out}{Gene/protein name of the connected protein}
#'   \item{score}{Interaction score of the selected connection}
#' }
#'
#' @export

place_singletons <- function( df_out, interactions, mapped_proteins, min_score ){

  ens_singletons <- df_out[df_out$cluster=='not_in_cluster', ]$string_name
  ens_in_cluster <- df_out[!df_out$cluster=='not_in_cluster', ]$string_name


  results <- list()
  for( ens_i in ens_singletons ){

    ind_1 <- interactions$from == ens_i | interactions$to == ens_i
    ind_2 <- interactions$from %in% ens_in_cluster | interactions$to %in% ens_in_cluster


    interactions_filtered <- interactions[ind_1 & ind_2, ]

    # NEW: apply score threshold
    interactions_filtered <- interactions_filtered[
      !is.na(interactions_filtered$score) &
        interactions_filtered$score >= min_score,
    ]


    if( nrow(interactions_filtered) == 0 ){
      next
    }


    ind_max <- which.max(interactions_filtered$score)


    score <- ifelse(length(ind_max)>0, interactions_filtered$score[ind_max], NA)
    connection_to <- ifelse(length(ind_max)>0, interactions_filtered$to[ind_max], NA)
    connection_from <- ifelse(length(ind_max)>0, interactions_filtered$from[ind_max], NA)

    connection <- ifelse(connection_to == ens_i, connection_from, connection_to)
    con_gn <- mapped_proteins[mapped_proteins$STRING_id==connection, ]$protein
    connection_gene_out <- ifelse(length(con_gn)>0, con_gn, ens_i)


    gn <- mapped_proteins[mapped_proteins$STRING_id==ens_i, ]$protein
    gene_out<- ifelse(length(gn)>0, gn, ens_i)

    results[[ length(results) + 1 ]] <- data.frame(ens_id = ens_i, gene_out, connection, connection_gene_out, score)

  }

  df_singleton <- if(length(results) > 0){

  do.call(rbind, results)

} else {

  data.frame(
    ens_id = character(0),
    gene_out = character(0),
    connection = character(0),
    connection_gene_out = character(0),
    score = numeric(0)
  )

}

return(na.omit(df_singleton))
}

