## Longitudinal pathway network for NICU-NEC


# Setup
setwd("~/Documents")
library(MicrobiomeDB, quietly = TRUE)
library(igraph, quietly = TRUE)
library(ggplot2)
library(stringr)

# Get the daily baby genus data
species_collection <- getCollection(
  microbiomeData::NICU_NEC,
  "Shotgun metagenomics Genus (Relative taxonomic abundance analysis)",
  continuousMetadataOnly = FALSE
)
pathway_collection <- getCollection(
  microbiomeData::NICU_NEC,
  "Shotgun metagenomics Metagenome enzyme pathway abundance data",
  continuousMetadataOnly = FALSE
)

# Subset by diagnosis day

# We need to filter by ourselves, which means deconstructing the collection
sampleMetadata <- getSampleMetadata(species_collection)
speciesAssayData <- microbiomeComputations::getAbundances(species_collection)
pathwayAssayData <- microbiomeComputations::getAbundances(pathway_collection)
ancestorIdColNames <- species_collection@ancestorIdColumns # Remove me from calculations
recordIColName <- species_collection@recordIdColumn # Use me to match data

diagnosis_day <- sort(unique(sampleMetadata$days_of_period_nec_diagnosed_days))
diagnosis_day <- -14:1


# Create a list of correlation graphs, one for each age in ages.
graph_list <- lapply(diagnosis_day, function(day) {
  
  # Subset to the appropriate abundances for this age
  day_samples <- sampleMetadata$Sample_Id[!is.na(sampleMetadata$days_of_period_nec_diagnosed_days) & sampleMetadata$days_of_period_nec_diagnosed_days == day]
  if (length(day_samples) < 5) {print(day); return(list())}
  day_species_abundances <- speciesAssayData[which(speciesAssayData$Sample_Id %in% day_samples), ]
  day_pathways_abundances <- pathwayAssayData[which(pathwayAssayData$Sample_Id %in% day_samples), ]
  print(day)
  
  # Make a new collection so we can send it through the correlation pipeline
  day_species_collection <- microbiomeComputations::AbundanceData(
    data=day_species_abundances,
    name=paste0("day_", day),
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames)
  
  day_pathways_collection <- microbiomeComputations::AbundanceData(
    data=day_pathways_abundances,
    name=paste0("day_", day),
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames)
  
  pathway_vs_species <- correlation(day_species_collection, day_pathways_collection, method = 'spearman')
  
  
  # Filter and extract data.table. IMPORTANT, this filtering here is our definition
  # of "correlated". To create the network, we must turn correlation into a binary
  # thing - either a taxon is correlated with a pathway or it isn't. To make this
  # more clear, when i mean correlated and passing filters, i'll write "correlated".
  corr_stats <- data.table::setDT(pathway_vs_species@statistics@statistics)
  corr_stats_filtered <- corr_stats[corr_stats$correlationCoef >= 0.5 & corr_stats$pValue <= 0.05, ]
  
  
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
  network <- igraph::graph_from_data_frame(network_df, directed=FALSE)
  # We may have many duplicate edges. If two nodes have k edges between them, it
  # means they are "correlated" with the same k pathways. 
  # We are not using that information at this time, so we'll simplify the graph
  # by removing duplicate edges.
  simplified_network <- simplify(
    network,
    remove.multiple = TRUE,
    remove.loops = FALSE
  )
  
  E(simplified_network)$sharedPathway <- network_df$sharedPathway
  
  V(simplified_network)$color <- unlist(lapply(V(simplified_network)$name, function(v) {ifelse(v %in% c("Clostridium", "Klebsiella"), "blue", "black")}))
  V(simplified_network)$label.color <- V(simplified_network)$color
  
  igraph::plot.igraph(
    simplified_network,
    arrow.mode=0,
    vertex.label.dist=3,
    vertex.label.degree=0,
    vertex.size=2,
    main=paste("Pathway network day", day)
    # edge.width=E(simplified_network)$weight/10
  )
  return(simplified_network)
})

## Now can do network stats on the graph list
cluster_results <- data.frame(nclusters=numeric(), max_size=numeric())
for (g in graph_list){
  clustering <- cluster_infomap(g)
  cluster_results <- rbind(cluster_results, data.frame(
    nclusters=length(clustering),
    max_size=max(unlist(lapply(seq_along(clustering),function(x){length(clustering[[x]])})))
  ))
}

n_pathways <- unlist(lapply(graph_list, function(g) {
  included_pathways <- unique(unlist(strsplit(unique(E(g)$sharedPathway), ",, ")))
  return(length(included_pathways))
}))


gg_data <- data.frame(diagnosis_day, cluster_results$nclusters, cluster_results$max_size, n_pathways)

# Plot
ggplot(gg_data, aes(x=diagnosis_day, y=cluster_results.nclusters)) +
  geom_line()

ggplot(gg_data, aes(x=diagnosis_day, y=cluster_results.max_size)) +
  geom_line()

ggplot(gg_data, aes(x=diagnosis_day, y=n_pathways)) +
  geom_line()




