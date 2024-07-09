## Longitudinal networks
library(ggplot2)

# Daily baby
DB_genus <- getCollection(
  microbiomeData::DailyBaby,
  "16S (V4) Genus (Relative taxonomic abundance analysis)",
  continuousMetadataOnly = FALSE
)

sampleMetadata <- getSampleMetadata(DB_genus)
assayData <- microbiomeComputations::getAbundances(DB_genus)
ancestorIdColNames <- DB_genus@ancestorIdColumns # Remove me from calculations
recordIColName <- DB_genus@recordIdColumn # Use me to match data

# Subset by age

ages = 1:10

graph_list <- lapply(ages, function(x) {
  
  age_samples <- sampleMetadata$Sample_Id[sampleMetadata$age_days == x]
  age_abundances <- assayData[which(assayData$Sample_Id %in% age_samples), ]
  
  age_collection <- microbiomeComputations::AbundanceData(
    data=age_abundances,
    name=toString(x),
    recordIdColumn = recordIColName,
    ancestorIdColumns = ancestorIdColNames)
  
  
  age_correlation <- selfCorrelation(age_collection, method = 'pearson')
  corr_stats <- data.table::setDT(age_correlation@statistics@statistics)
  filtered_corr_stats <- corr_stats[corr_stats$correlationCoef >= 0.3 & corr_stats$pValue <= 0.05, ]
  age_graph <- igraph::graph_from_data_frame(filtered_corr_stats, directed=FALSE)

  igraph::plot.igraph(
    age_graph,
    arrow.mode=0,
    vertex.color="white",
    vertex.label.dist=1,
    vertex.label.color="black",
    vertex.label.degree=0,
    vertex.size=2,
    main=paste("day", x, sep=' ')
  )
  return(age_graph)
})

## Now can do network stats on the graph list
# NOTE the graph is currenty unweighted so this is just avg degree
avg_strength <- unlist(lapply(graph_list, function(g) {
  mean(strength(g))
  # transitivity(g)
}))


gg_data <- data.frame(ages,avg_strength)

# Plot
ggplot(gg_data, aes(x=ages, y=avg_strength)) +
  geom_line()
