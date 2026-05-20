#' Add singleton proteins to the interaction graph
#'
#' Extends a filtered protein interaction network by adding singleton proteins
#' and connecting them to their strongest interaction partners. Singleton edges
#' are added with distinct styling to differentiate them from original network edges.
#'
#' The function also computes a layout for visualization and annotates nodes
#' with gene names.
#'
#' This corresponds to the network augmentation and visualization preparation
#' step in the ProteoNet pipeline.
#'
#' @param interaction_graph_filtered An igraph object representing the filtered network
#' @param df_singleton Data frame of singleton assignments (e.g. from \code{place_singletons})
#' @param df_out Data frame of subnetwork membership assignments
#' @param mapped_proteins Data frame mapping STRING IDs to gene/protein names
#' @param protein_highlight Character vector of proteins to highlight (currently unused)
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{graph}{An igraph object including singleton nodes and edges}
#'   \item{layout}{A ggraph layout object with node positions and annotations}
#' }
#'
#' @export


add_singletons_to_graph <- function(interaction_graph_filtered, df_singleton, df_out, mapped_proteins, protein_highlight){

  df_out <- df_out[!is.na(df_out$gene), ]

  if(nrow(df_singleton) == 0){

    layout_fr <- ggraph::create_layout(
      interaction_graph_filtered,
      layout = "fr"
    )

    gene_name <- sapply(
      layout_fr$name,
      function(str_i){

        mapped_proteins$protein[
          mapped_proteins$STRING_id == str_i
        ][1]

      }
    )

    layout_fr$gene_name <- gene_name

    return(list(
      layout = layout_fr,
      graph = interaction_graph_filtered
    ))

  }



  new_vertices <- c()
  new_community <- c()

  for(ens_i in df_singleton$ens_id){

    new_vertices <- c(new_vertices, ens_i)

    df_i <- df_singleton[df_singleton$ens_id == ens_i, ]

    new_community <- c(new_community, as.numeric(df_out[df_out$string_name == df_i$connection, "cluster"]))

  }

  interaction_graph_with_singletons <- igraph::add_vertices(
    interaction_graph_filtered,
    nv   = length(new_vertices),
    name = new_vertices,
    community = new_community
  )


  for(i in seq_len(nrow(df_singleton))){

    interaction_graph_with_singletons <- igraph::add_edges(
      interaction_graph_with_singletons,
      c(df_singleton$ens_id[i], df_singleton$connection[i]),
      attr = list(alternative_score = df_singleton$score[i],
                  linestyle  = "dashed" )
    )

  }


  igraph::E(interaction_graph_with_singletons)$linestyle[
    is.na(igraph::E(interaction_graph_with_singletons)$linestyle)
  ] <- "solid"

  layout_fr <- ggraph::create_layout(
    interaction_graph_with_singletons,
    layout = "fr"
  )


  gene_name <- c()
  for(str_i in layout_fr$name){
    gene_name <- c(gene_name, mapped_proteins$protein[mapped_proteins$STRING_id == str_i][1])
  }

  layout_fr$gene_name <- gene_name



  return(list(layout = layout_fr, graph = interaction_graph_with_singletons))
}
