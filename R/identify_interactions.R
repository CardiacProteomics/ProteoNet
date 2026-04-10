#' Build a protein-protein interaction network from STRING
#'
#' Constructs a protein interaction network using a locally downloaded
#' STRING database for the specified species.
#'
#' This corresponds to the network construction step in the ProteoNet pipeline.
#'
#' @param proteins Character vector of gene or protein identifiers
#' @param string_database_location Path to the folder containing STRING data
#' @param species Numeric species identifier (e.g. 9606 for Homo sapiens)
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{interactions}{Object with protein-protein interactions}
#'   \item{string_object}{STRING database object used for mapping}
#'   \item{mapping}{Data frame mapping STRING IDs to gene symbols}
#' }
#' @export


identify_interactions <- function( proteins, string_database_location, species ){
  
  string_db <- STRINGdb::STRINGdb$new(version="12", species=species, score_threshold=0, input_directory=string_database_location, link_data='detailed')
  mapped_proteins <- string_db$map(data.frame(protein = proteins), "protein", removeUnmappedRows = TRUE)
  all_interactions_in <- string_db$get_interactions(mapped_proteins$STRING_id)
  interactions <- all_interactions_in |> dplyr::distinct(from, to, .keep_all = TRUE)
  
  interactions$score <- 1 - (1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$database/1000)
  
  interactions$alternative_score <- 1 -(1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$neighborhood/1000)*
    (1-interactions$fusion/1000)*(1-interactions$cooccurence/1000)*(1-interactions$textmining/1000)
  
  map_temp <- merge(
    mapped_proteins,
    string_db$proteins,
    by.x = "STRING_id",
    by.y = "protein_external_id"
  )
  
  mapped_proteins <- map_temp |>
    dplyr::select(-preferred_name, -protein_size)
  
  list(interactions=interactions, string_db = string_db, mapped_proteins = mapped_proteins)
}

