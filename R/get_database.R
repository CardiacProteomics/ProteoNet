

##https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp


get_database <- function( database, folder ){
  
  if(database == 'gobp'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/c5.go.bp.v2026.1.Hs.symbols.gmt"))
  }
  
  if(database == 'gocc'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/c5.go.cc.v2026.1.Hs.symbols.gmt"))
  }
  
  if(database == 'gomf'){
    pathways <- fgsea::gmtPathways(paste0(folder, "/Gene_set_databases/m5.go.mf.v2024.1.Mm.symbols.gmt"))
  }
  
  if(database == 'reactome'){
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

