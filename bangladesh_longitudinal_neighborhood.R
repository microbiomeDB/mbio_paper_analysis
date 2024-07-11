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



# Create a list of correlation graphs, one for each age in ages. The original paper
# separated them by month, so let's do the same
ages <- floor(min(sampleMetadata$age_months)):ceiling(max(sampleMetadata$age_months))

graph_list <- lapply(1:(length(ages)-1), function(x) {
  
  # Subset to the appropriate abundances for this age
  age_samples <- sampleMetadata$Sample_Id[(sampleMetadata$age_months >= ages[x] & sampleMetadata$age_months < ages[x+1])]

  # Bail if there aren't enough samples
  if (length(age_samples) < min_samples_required) {
    return(c())
  } else {

    age_abundances <- assayData[which(assayData$Sample_Id %in% age_samples), ]
    
    # Make a new collection so we can send it through the correlation pipeline
    age_collection <- microbiomeComputations::AbundanceData(
      data=age_abundances,
      name=paste("Age", ages[x], "to", ages[x+1], "months"),
      recordIdColumn = recordIColName,
      ancestorIdColumns = ancestorIdColNames)
    
    age_correlation <- selfCorrelation(age_collection, method = 'pearson')
  
    # Extract results. We want to do this ourselves to have more control over
    # how we create the igraph object.
    corr_stats <- data.table::setDT(age_correlation@statistics@statistics)
    filtered_corr_stats <- corr_stats[abs(corr_stats$correlationCoef) >= 0.3 & corr_stats$pValue <= 0.05, ]
    age_graph <- igraph::graph_from_data_frame(filtered_corr_stats, directed=FALSE)
  
    # Now we want just the neighborhood of the node we care about
    # Note an "ego graph" means "node neighborhood". Order=1 means only keep 
    # the direct neighbors of the node of interest. Returns a list of graphs
    cool_taxon_neighborhood <- make_ego_graph(age_graph, order=1, nodes=which(V(age_graph)$name == cool_taxon))
    V(cool_taxon_neighborhood[[1]])$color <- "white"
    V(cool_taxon_neighborhood[[1]])$color[V(cool_taxon_neighborhood[[1]])$name == cool_taxon] <- "red"
    
    # Plot! But only if there are any vertices
    if (length(cool_taxon_neighborhood)>0 && length(V(cool_taxon_neighborhood[[1]])) > 0) {
      igraph::plot.igraph(
        cool_taxon_neighborhood[[1]],
        arrow.mode=0,
        # vertex.color="white",
        vertex.label.dist=10,
        vertex.label.color="black",
        vertex.label.degree=0,
        vertex.size=10,
        main=paste("Age", ages[x], "to", ages[x+1], "months")
      )

    }
    age_graph$age <- ages[x]
    
    return(age_graph)
  }
})

## Now can do network stats on the graph list
# Calculate cool_taxon's degree for each graph
# Note that the plot looks weird because cool_taxon gets degree 0 if there
# weren't enough samples, or if it didn't actually have any correlation
# friends that passed the thresholds
cool_degrees <- unlist(lapply(graph_list[1:60], function(g) {
  if (length(g) > 0) {
    node_degrees <- degree(g)
    cool_degree <- node_degrees[which(names(node_degrees) == cool_taxon)]
    if (length(cool_degree) ==0) {
      return(0)
    }
    return(cool_degree)
  } else {
    return(0)
  }
}))


gg_data <- data.frame(ages=unlist(lapply(graph_list[1:60], function(x) x$age)),cool_degrees)

# Plot
ggplot(gg_data, aes(x=ages, y=cool_degrees)) +
  geom_line() +
  geom_point() +
  labs(x="Age (months)", y=paste("Degree of", cool_taxon), title=paste("Degree of", cool_taxon, "across early life")) +
  theme_bw()


