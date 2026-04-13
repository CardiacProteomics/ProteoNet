#' Load gene set database for enrichment analysis
#'
#' Loads a gene set database in GMT format from a local folder using a
#' predefined database identifier (e.g. GO, Reactome).
#'
#' The function maps shorthand database names to specific GMT files and
#' returns pathway definitions suitable for over-representation analysis.
#'
#' @param database Character string specifying the database to load.
#'   Supported options include:
#'   \code{"gobp"}, \code{"gocc"}, \code{"gomf"},
#'   \code{"reactome"}, \code{"reactome_human"},
#'   \code{"gobp_human"}, \code{"gocc_human"}
#' @param folder Path to the base folder containing the
#'   \code{Gene_set_databases} directory with GMT files
#'
#' @return A list of pathways, where each element is a character vector
#'   of genes belonging to a gene set (as returned by \code{fgsea::gmtPathways})
#'
#' @export


get_database <- function( database, folder ){

  if(database == 'gobp_mouse'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/m5.go.bp.v2024.1.Mm.symbols.gmt"))
  }

  if(database == 'gocc_mouse'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/m5.go.cc.v2024.1.Mm.symbols.gmt"))
  }

  if(database == 'gomf_mouse'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/m5.go.mf.v2024.1.Mm.symbols.gmt"))
  }

  if(database == 'reactome_mouse'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/m2.cp.reactome.v2024.1.Mm.symbols.gmt"))
  }

  if(database == 'reactome_human'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
  }

  if(database == 'gobp_human'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/c5.go.bp.v2025.1.Hs.symbols.gmt"))
  }

  if(database == 'gocc_human'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/c5.go.cc.v2025.1.Hs.symbols.gmt"))
  }

  return(pathways)
}

