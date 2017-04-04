loadGeneticMap <- function(chrom,filename){
  genetic.map.db <- dbConnect( dbDriver( "SQLite" ), filename )
  gm = dbGetQuery(
    genetic.map.db,
    sprintf(
      "SELECT * FROM GeneticMap WHERE chromosome == '%s' AND position >= %d AND position <= %d ORDER BY position",
      chrom, range(snp.pos)[1], range(snp.pos)[2]
    )
  )
  
  
  gm$COMBINED_LD_rate = NA
  gm$COMBINED_LD_rate[-1] = ( gm$COMBINED_LD[-1] - gm$COMBINED_LD[-nrow(gm)] ) / ( gm$position[-1] - gm$position[-nrow(gm)] )
  gm$YRI_LD_rate = NA
  gm$YRI_LD_rate[-1] = ( gm$YRI_LD[-1] - gm$YRI_LD[-nrow(gm)] ) / ( gm$position[-1] - gm$position[-nrow(gm)] )
  gm$African_Enriched_rate = NA
  gm$African_Enriched_rate[-1] = ( gm$African_Enriched[-1] - gm$African_Enriched[-nrow(gm)] ) / ( gm$position[-1] - gm$position[-nrow(gm)] )
  return(gm)
}
