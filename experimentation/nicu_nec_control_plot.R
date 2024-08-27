## NICU-NEC Control plot

setwd("~/Documents")
library(MicrobiomeDB, quietly = TRUE)
library(igraph, quietly = TRUE)
library(ggplot2)
library(stringr)
source("nicu_nec_analysis/utils.R")

# Retrieve data
species_collection <- getCollection(
  microbiomeData::NICU_NEC,
  "Shotgun metagenomics Species (Relative taxonomic abundance analysis)",
  continuousMetadataOnly = FALSE
)
pathway_collection <- getCollection(
  microbiomeData::NICU_NEC,
  "Shotgun metagenomics Metagenome enzyme pathway abundance data",
  continuousMetadataOnly = FALSE
)

# We need to filter by ourselves, which means deconstructing the collection
sampleMetadata <- getSampleMetadata(species_collection)
speciesAssayData <- microbiomeComputations::getAbundances(species_collection)
pathwayAssayData <- microbiomeComputations::getAbundances(pathway_collection)
ancestorIdColNames <- species_collection@ancestorIdColumns # Remove me from calculations
recordIColName <- species_collection@recordIdColumn # Use me to match data

# We want to create a plot of network stats by age for the control patients. 
ages <- 0:70

graph_list <- lapply(ages, function(age) {
  
  # So, we need to first filter to only the control patients.
  # Making the assumption that samples with NA for diagnosis day are control patients.
  age_samples <- sampleMetadata$Sample_Id[is.na(sampleMetadata$days_of_period_nec_diagnosed_days) & sampleMetadata$age_days == age]
  
  if (length(age_samples) < 5) {print(paste("Not enough samples for age", age)); return(list())}
  
  age_species_abundances <- speciesAssayData[which(speciesAssayData$Sample_Id %in% age_samples), ]
  age_pathways_abundances <- pathwayAssayData[which(pathwayAssayData$Sample_Id %in% age_samples), ]
  
  # Make a new collection so we can send it through the correlation pipeline
  age_species_collection <- microbiomeComputations::AbundanceData(
    data=age_species_abundances,
    name=paste0("age_", age),
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames)

  age_pathways_collection <- microbiomeComputations::AbundanceData(
    data=age_pathways_abundances,
    name=paste0("age_", age),
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames
  )

  shared_pathway_network <- createSharedPathwayNetwork(age_species_collection, age_pathways_collection, method="spearman", corrCoeffThreshold = 0.5, pValueThreshold = 0.05)
  
  # igraph::plot.igraph(
  #   shared_pathway_network,
  #   arrow.mode=0,
  #   vertex.label.dist=3,
  #   vertex.label.degree=0,
  #   vertex.size=2,
  #   main=paste("Conrol pathway network age", age)
  #   # edge.width=E(simplified_network)$weight/10
  # )

  return(shared_pathway_network)

})


## Plot!
## Now can do network stats on the graph list
cluster_results <- data.frame(nclusters=numeric(), max_size=numeric())
for (g in graph_list){
  if(length(g) == 0) {
    cluster_results <- rbind(cluster_results, data.frame(
      nclusters=0,
      max_size=0
    ))
    print(nrow(cluster_results))
  } else {
    clustering <- cluster_infomap(g)
    cluster_results <- rbind(cluster_results, data.frame(
      nclusters=length(clustering),
      max_size=max(unlist(lapply(seq_along(clustering),function(x){length(clustering[[x]])})))
    ))
  }
}

n_pathways <- unlist(lapply(graph_list, function(g) {
  if(length(g) == 0) return(0)
  included_pathways <- unique(unlist(strsplit(unique(E(g)$sharedPathway), ",, ")))
  return(length(included_pathways))
}))


gg_data <- data.frame(ages, cluster_results$nclusters, cluster_results$max_size, n_pathways)

# Plot
ggplot(gg_data, aes(x=ages, y=cluster_results.nclusters)) +
  geom_line()

ggplot(gg_data, aes(x=ages, y=cluster_results.max_size)) +
  geom_line()

ggplot(gg_data, aes(x=ages, y=n_pathways)) +
  geom_line()


## Can we see anything different with the nec patients? Unlikely but let's try
# So maybe it was a bit different, but it's pretty inconclusive.

