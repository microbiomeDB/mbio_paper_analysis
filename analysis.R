## R analysis for mbio

# Setup
setwd("~/Documents")
library(MicrobiomeDB, quietly = TRUE)
library(igraph, quietly = TRUE)

# Run microbiomeData::getCuratedDatasetNames() to get all names of studies
# or getCollectionNames(microbiomeData::HMP_V1V3) to get all collection names
HMP_V1V3_genus <- getCollection(
  microbiomeData::HMP_V1V3,
  "16S (V1-V3) Genus (Relative taxonomic abundance analysis)",
  continuousMetadataOnly = FALSE
  )

sampleMetadata <- getSampleMetadata(HMP_V1V3_genus)
assayData <- microbiomeComputations::getAbundances(HMP_V1V3_genus)
ancestorIdColNames <- HMP_V1V3_genus@ancestorIdColumns # Remove me from calculations
recordIColName <- HMP_V1V3_genus@recordIdColumn # Use me to match data

# I want to make four subsets based on host body habitat
habitats <- unique(sampleMetadata$host_body_habitat)

graph_list <- lapply(habitats, function(x) {
  
  habitat_samples <- sampleMetadata$Sample_Id[sampleMetadata$host_body_habitat == x]
  habitat_abundances <- assayData[which(assayData$Sample_Id %in% habitat_samples), ]
  
  habitat_collection <- microbiomeComputations::AbundanceData(
    data=habitat_abundances,
    name=x,
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames)
  
  
  habitat_correlation <- selfCorrelation(habitat_collection, method = 'pearson')
  corr_stats <- data.table::setDT(habitat_correlation@statistics@statistics)
  filtered_corr_stats <- corr_stats[corr_stats$correlationCoef >= 0.3 & corr_stats$pValue <= 0.05, ]
  habitat_graph <- igraph::graph_from_data_frame(filtered_corr_stats, directed=FALSE)

  igraph::plot.igraph(
    habitat_graph,
  arrow.mode=0,
  vertex.color="white",
  vertex.label.dist=1,
  vertex.label.color="black",
  vertex.label.degree=0,
  vertex.size=2,
  main=x
  )
  return(habitat_graph)
})

