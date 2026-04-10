
library(igraph)
library(ggraph)

identify_interactions <- function( my_geneset, string_database_location, species ){
  
  string_db <- STRINGdb::STRINGdb$new(version="12", species=species, score_threshold=0, input_directory=string_database_location, link_data='detailed')
  mapped_proteins <- string_db$map(data.frame(protein = my_geneset), "protein", removeUnmappedRows = TRUE)
  all_interactions_in <- string_db$get_interactions(mapped_proteins$STRING_id)
  interactions <- all_interactions_in %>% dplyr::distinct(from, to, .keep_all = TRUE)
  
  interactions$score <- 1 - (1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$database/1000)
  
  interactions$alternative_score <- 1 -(1-interactions$coexpression/1000)*(1 - interactions$experimental/1000)*(1 - interactions$neighborhood/1000)*
    (1-interactions$fusion/1000)*(1-interactions$cooccurence/1000)*(1-interactions$textmining/1000)
  
  map_temp <- merge(
    mapped_proteins,
    string_db$proteins,
    by.x = "STRING_id",
    by.y = "protein_external_id"
  )
  
  mapped_proteins <- map_temp %>%
    dplyr::select(-preferred_name, -protein_size)
  
  list(interactions=interactions, string_db = string_db, mapped_proteins = mapped_proteins)
}

