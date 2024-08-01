## R analysis for mbio

# Setup
setwd("~/Documents")
library(MicrobiomeDB, quietly = TRUE)
library(igraph, quietly = TRUE)

# Run microbiomeData::getCuratedDatasetNames() to get all names of studies
# or getCollectionNames(microbiomeData::HMP_V1V3) to get all collection names
MALED_genus <- getCollection(
  microbiomeData::NICU_NEC,
  "Shotgun metagenomics Species (Relative taxonomic abundance analysis)",
  continuousMetadataOnly = FALSE
)

sampleMetadata <- getSampleMetadata(MALED_genus)
assayData <- microbiomeComputations::getAbundances(MALED_genus)
ancestorIdColNames <- MALED_genus@ancestorIdColumns # Remove me from calculations
recordIColName <- MALED_genus@recordIdColumn # Use me to match data

# I want to make subsets based on country
countries <- unique(sampleMetadata$amnionicity)
countries <- countries[2:4]


graph_list <- lapply(countries, function(x) {
  
  case_samples <- sampleMetadata$Sample_Id[sampleMetadata$amnionicity == x]
  case_abundances <- assayData[which(assayData$Sample_Id %in% case_samples), ]
  
  case_collection <- microbiomeComputations::AbundanceData(
    data=case_abundances,
    name=x,
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames)
  
  
  case_correlation <- selfCorrelation(case_collection, method = 'sparcc')
  corr_stats <- data.table::setDT(case_correlation@statistics@statistics)
  filtered_corr_stats <- corr_stats[abs(corr_stats$correlationCoef) >= 0.25 & corr_stats$pValue <= 0.05, ]
  case_graph <- igraph::graph_from_data_frame(filtered_corr_stats, directed=FALSE)

  igraph::plot.igraph(
    case_graph,
  arrow.mode=0,
  vertex.color="white",
  vertex.label.dist=1,
  vertex.label.color="black",
  vertex.label.degree=0,
  vertex.size=degree(case_graph)*4,
  main=x
  )
  return(case_graph)
})


# Are any nodes overlapping?
brazil_nodes <- V(graph_list[[1]])$name
india_nodes <- V(graph_list[[2]])$name
both_nodes <- intersect(brazil_nodes, india_nodes)

# There's a good intersection here. Make into single weighted network to create
# a combined layout. Then use this layout to plot and highlight nodes in 
# both networks
brazil_df <- igraph::as_long_data_frame(graph_list[[1]])
india_df <- igraph::as_long_data_frame(graph_list[[2]])
combined_df <- rbind(brazil_df, india_df)
# Simplify combined_df to have one edge per pair of nodes and add a column for number of times seen
combined_df$occurrences <- ave(combined_df$from_name, combined_df$from_name, combined_df$to_name, FUN = length)
# igraph is very particular about how the input df is set up, so let's redo it
combined_df <- combined_df[, c("from_name", "to_name", "occurrences")]

combined_graph <- igraph::graph_from_data_frame(combined_df, directed=FALSE)
E(combined_graph)$weight <- combined_df$occurrences
brazil_graph <- graph_list[[1]]
india_graph <- graph_list[[2]]
coords <- as.data.frame(layout_nicely(combined_graph))
coords$node_name <- V(combined_graph)$name


# Plot brazil
V(brazil_graph)$color <- unlist(lapply(V(brazil_graph)$name, function(v) {ifelse(v %in% both_nodes, "blue", "black")}))
V(brazil_graph)$label.color <- V(brazil_graph)$color
brazil_coords <- as.matrix(coords[match(V(brazil_graph)$name, coords$node_name), 1:2])
igraph::plot.igraph(
  brazil_graph,
  arrow.mode=0,
  vertex.label.dist=1,
  vertex.label.degree=0,
  vertex.size=10,
  main="brazil",
  layout=brazil_coords,
  rescale=F,xlim=c(-6,3),ylim=c(-8,3)
)



V(india_graph)$color <- unlist(lapply(V(india_graph)$name, function(v) {ifelse(v %in% both_nodes, "blue", "black")}))
V(india_graph)$label.color <- V(india_graph)$color
india_coords <- as.matrix(coords[match(V(india_graph)$name, coords$node_name), 1:2])
igraph::plot.igraph(
  india_graph,
  arrow.mode=0,
  vertex.label.dist=1,
  vertex.label.degree=0,
  vertex.size=10,
  main='india',
  layout=india_coords,
  rescale=F,xlim=c(-6,3),ylim=c(-8,3)
)

