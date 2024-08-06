## Function to return the shared pathway network

createSharedPathwayNetwork <- function(taxaCollection, pathwayCollection, method, corrCoeffThreshold, pValueThreshold) {
  
  # Compute the correlation
  pathway_vs_species <- correlation(taxaCollection, pathwayCollection, method = method)
  
  
  # Filter and extract data.table. IMPORTANT, this filtering here is our definition
  # of "correlated". To create the network, we must turn correlation into a binary
  # thing - either a taxon is correlated with a pathway or it isn't. To make this
  # more clear, when i mean correlated and passing filters, i'll write "correlated".
  corr_stats <- data.table::setDT(pathway_vs_species@statistics@statistics)
  # Assume a negative threshold means <= and a positive one means >=
  if (corrCoeffThreshold < 0) {
    corr_stats_filtered <- corr_stats[corr_stats$correlationCoef <= corrCoeffThreshold & corr_stats$pValue <= pValueThreshold, ]
  } else {
    corr_stats_filtered <- corr_stats[corr_stats$correlationCoef >= corrCoeffThreshold & corr_stats$pValue <= pValueThreshold, ]
  }
  
  
  # Now what we need to do is turn this list of "correlations" into a network
  # I'll use for loops for clarity.
  # Initialize empty data.frame
  network_df <- data.frame(
    taxonA = character(),
    taxonB = character(),
    sharedPathway = character(),
    stringsAsFactors = FALSE,
    weight = numeric()
  )
  
  pathways <- as.character(unique(corr_stats_filtered$data2))
  for (pathway in pathways) {
    taxa_correlated_with_pathway <- corr_stats_filtered$data1[which(corr_stats_filtered$data2 == pathway)]
    # We need at least two taxa to be correlated with the pathway to make any edges!
    if (length(taxa_correlated_with_pathway) > 1){
      # Add all pairs of taxa to the network data frame
      for (i in 1:(length(taxa_correlated_with_pathway)-1)) {
        taxon <- taxa_correlated_with_pathway[i]
        for (j in (i+1):length(taxa_correlated_with_pathway)) {
          # Does this edge already exist?
          if(nrow(network_df[network_df$taxonA == taxon & network_df$taxonB == taxa_correlated_with_pathway[j], ]) > 0) {
            # If so, increment the weight
            network_df[network_df$taxonA == taxon & network_df$taxonB == taxa_correlated_with_pathway[j], ]$weight <- network_df[network_df$taxonA == taxon & network_df$taxonB == taxa_correlated_with_pathway[j], ]$weight + 1
            # Append this pathway to the pathway list
            network_df[network_df$taxonA == taxon & network_df$taxonB == taxa_correlated_with_pathway[j], ]$sharedPathway <- paste(network_df[network_df$taxonA == taxon & network_df$taxonB == taxa_correlated_with_pathway[j], ]$sharedPathway, pathway, sep=",, ")
          } else {
            # If not, add the edge
            network_df <- rbind(network_df, data.frame(taxonA = taxon, taxonB = taxa_correlated_with_pathway[j], sharedPathway=pathway, weight=1, stringsAsFactors = FALSE))
          }
        }
      }
    } else if (length(taxa_correlated_with_pathway) == 1) {
      # Then there's just a lonely node. Let's add it with a self-loop to keep track of 
      # all microbes doing things in the system
      taxon <- taxa_correlated_with_pathway[1]
      # Does this edge already exist?
      if(nrow(network_df[network_df$taxonA == taxon & network_df$taxonB == taxon, ]) > 0) {
        # If so, increment the weight
        network_df[network_df$taxonA == taxon & network_df$taxonB == taxon, ]$weight <- network_df[network_df$taxonA == taxon & network_df$taxonB == taxon, ]$weight + 1
        # Append this pathway to the pathway list
        network_df[network_df$taxonA == taxon & network_df$taxonB == taxon, ]$sharedPathway <- paste(network_df[network_df$taxonA == taxon & network_df$taxonB == taxon, ]$sharedPathway, pathway, sep=",, ")
      } else {
        # If not, add the edge
        network_df <- rbind(network_df, data.frame(taxonA = taxon, taxonB = taxon, sharedPathway=pathway, weight=1, stringsAsFactors = FALSE))
      }
    }
  }
  
  # Now make an igraph network
  # Note, there's no need to simplify the resulting network because above we
  # ensured that there are no duplicate edges, and also we want to keep all the self-loops
  network <- igraph::graph_from_data_frame(network_df, directed=FALSE)

  
  E(network)$sharedPathway <- network_df$sharedPathway

  return(network)
}

## Rescale values to be between min_value and max_value
rescale <- function(x, a=1, b=20, min_value=0, max_value=1) {
  return((b - a) * (x - min_value) / (max_value - min_value) + a)
}