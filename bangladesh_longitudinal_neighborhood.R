## Longitudinal neighborhood analysis of a txon of interest
library(ggplot2)

# Get the bangladesh genus data
Bangladesh_genus <- getCollection(
  microbiomeData::Bangladesh,
  "16S (V4) Genus (Relative taxonomic abundance analysis)",
  continuousMetadataOnly = FALSE
)

# Subset by age
min_samples_required <- 10
cool_taxon <- 'Bifidobacterium'

# We need to filter by ourselves, which means deconstructing the collection
sampleMetadata <- getSampleMetadata(Bangladesh_genus)
assayData <- microbiomeComputations::getAbundances(Bangladesh_genus)
ancestorIdColNames <- Bangladesh_genus@ancestorIdColumns # Remove me from calculations
recordIColName <- Bangladesh_genus@recordIdColumn # Use me to match data



# Create a list of correlation graphs, one for each age in ages.
ages <- sort(unique(sampleMetadata$age_months))
graph_list <- lapply(ages, function(x) {
  
  # Subset to the appropriate abundances for this age
  age_samples <- sampleMetadata$Sample_Id[sampleMetadata$age_months == x]

  # Bail if there aren't enough samples
  if (length(age_samples) < min_samples_required) {
    return(c())
  } else {

    age_abundances <- assayData[which(assayData$Sample_Id %in% age_samples), ]
    
    # Make a new collection so we can send it through the correlation pipeline
    age_collection <- microbiomeComputations::AbundanceData(
      data=age_abundances,
      name=toString(x),
      recordIdColumn = recordIColName,
      ancestorIdColumns = ancestorIdColNames)
    
    age_correlation <- selfCorrelation(age_collection, method = 'pearson')
  
    # Extract results. We want to do this ourselves to have more control over
    # how we create the igraph object.
    corr_stats <- data.table::setDT(age_correlation@statistics@statistics)
    filtered_corr_stats <- corr_stats[abs(corr_stats$correlationCoef) >= 0.10 & corr_stats$pValue <= 0.2, ]
    age_graph <- igraph::graph_from_data_frame(filtered_corr_stats, directed=FALSE)
  
    # Now we want just the neighborhood of the node we care about
    # Note an "ego graph" means "node neighborhood". Order=1 means only keep 
    # the direct neighbors of the node of interest. Returns a list of graphs
    cool_taxon_neighborhood <- make_ego_graph(age_graph, order=1, nodes=which(V(age_graph)$name == cool_taxon))
    
    # Plot! But only if there are any vertices
    if (length(cool_taxon_neighborhood)>0 && length(V(cool_taxon_neighborhood[[1]])) > 0) {
      igraph::plot.igraph(
        cool_taxon_neighborhood[[1]],
        arrow.mode=0,
        vertex.color="white",
        vertex.label.dist=1,
        vertex.label.color="black",
        vertex.label.degree=0,
        vertex.size=2,
        main=paste("month", x, sep=' ')
      )

    }
    
    return(age_graph)
  }
})

## Now can do network stats on the graph list
# Calculate cool_taxon's degree for each graph
# Note that the plot looks weird because cool_taxon gets degree 0 if there
# weren't enough samples, or if it didn't actually have any correlation
# friends that passed the thresholds
cool_degrees <- unlist(lapply(graph_list, function(g) {
  if (length(g) > 0) {
    node_degrees <- degree(g)
    cool_degree <- node_degrees[which(names(node_degrees) == cool_taxon)]
    return(cool_degree)
  } else {
    return(0)
  }
}))


gg_data <- data.frame(ages,cool_degrees)

# Plot
ggplot(gg_data, aes(x=ages, y=cool_degrees)) +
  geom_line()
