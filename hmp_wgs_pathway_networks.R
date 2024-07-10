# Shared pathway network
# Here we'll create a network of taxa. Taxa will be connected if the
# two taxa are both correlated with at least 1 pathway.
# A taxon connected to many other taxa in this network may be correlated with 
# lots of pathways, or it may be correlated with only pathways that are very
# popular among other taxa. 
# Groups of taxa have some overlapping pathways. 
# In particular, the true representation is a simplicial complex, becuase if 
# n taxa all correlate with a pathway, then all n taxa will
# be connected to each other.
# We can consider picking out pathways we care about and then looking at
# the filtered network that use these pathways.
# It would also make a nice heatmap


library(ggplot2)
library(reshape2)
library(viridis)

# Grab collections
HMP_MGX_species <- getCollection(microbiomeData::HMP_MGX, "Shotgun metagenomics Species (Relative taxonomic abundance analysis)")
HMP_MGX_pathways <- getCollection(microbiomeData::HMP_MGX, "Shotgun metagenomics Metagenome enzyme pathway abundance data" )

# Correlation
pathway_vs_species <- correlation(HMP_MGX_species, HMP_MGX_pathways, method = 'spearman')

# Filter and extract data.table. IMPORTANT, this filtering here is our definition
# of "correlated". To create the network, we must turn correlation into a binary
# thing - either a taxon is correlated with a pathway or it isn't. To make this
# more clear, when i mean correlated and passing filters, i'll write "correlated".
corr_stats <- data.table::setDT(pathway_vs_species@statistics@statistics)
corr_stats_filtered <- corr_stats[abs(corr_stats$correlationCoef) >= 0.80 & corr_stats$pValue <= 0.01, ]


# Now what we need to do is turn this list of "correlations" into a network
# I'll use for loops for clarity.
# Initialize empty data.frame
network_df <- data.frame(
  taxonA = character(),
  taxonB = character(),
  sharedPathway = character(),
  stringsAsFactors = FALSE
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
        # Add edge from taxon to taxa_correlated_with_pathway[j]
        network_df <- rbind(network_df, data.frame(taxonA = taxon, taxonB = taxa_correlated_with_pathway[j], sharedPathway=pathway, stringsAsFactors = FALSE))
      }
    }
  }
}

# Let's do a check to make sure we have the correct number of edges
# The incidence of each pathway gives us the size of the cliques, then we can
# calculate the number of edges from the cliques. There may be duplicates, but 
# they would be the same duplicates as above
pathway_frequencies <- as.data.frame(table(corr_stats_filtered$data2))$Freq
num_edges <- 0
for (pathway_frequency in pathway_frequencies) {
  # Calculate the number of edges in the clique
  num_edges <- num_edges + (pathway_frequency * (pathway_frequency - 1) / 2)
}
if (num_edges == nrow(network_df)) {
  print("Check complete")
} else {
  print("Hey there's probably a bug....")
}

# Now make an igraph network
network <- igraph::graph_from_data_frame(network_df, directed=FALSE)
# We may have many duplicate edges. If two nodes have k edges between them, it
# means they are "correlated" with the same k pathways. 
# We are not using that information at this time, so we'll simplify the graph
# by removing duplicate edges.
# simplified_network <- simplify(
#   network,
#   remove.multiple = TRUE
# )
E(network)$color <- as.factor(network_df$sharedPathway)

igraph::plot.igraph(
  g,
  arrow.mode=0,
  vertex.color="white",
  vertex.label.dist=1,
  vertex.label.color="black",
  vertex.label.degree=0,
  vertex.size=2,
  main="Pathway network"
)


## Want to know which pathways are connecting nodes?
corr_incidence <- reshape(corr_stats_filtered[, c("data1", "data2", "correlationCoef")], direction="wide", idvar="data1", timevar="data2")
# Quick n dirty heatmap
data <- melt(corr_incidence)
ggplot(data, aes(x = data1, y = variable, fill = value)) +
  geom_tile() +
  labs(title = "Correlation Heatmap. Abs(corr) >=0.8, pvalue <=0.01",
       x = "Taxon",
       y = "Pathway") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## Let's draw some nicer graphs
component_membership <- components(network)$membership
components_node_list <- lapply(seq(1:max(component_membership)), function(i) {
  return(V(network)$name[which(component_membership==i)])
})
connected_components_graphs <- lapply(seq(1:max(component_membership)), function(i) {
  g <- induced_subgraph(network, components_node_list[[i]], "create_from_scratch")
})


# Plot graph and legend
for(i in seq(1:max(component_membership))) {
  component_i <- connected_components_graphs[[i]]
  # Find the unique pathways
  component_i_pathways <- unique(E(component_i)$color)
  # Create a legend data.frame that maps the pathway to a color
  legend_data <- data.frame(
    pathway = component_i_pathways,
    color = viridis(length(component_i_pathways))
  )
  # Remap the edge colors to the legend colors
  E(component_i)$color <- legend_data$color[match(E(component_i)$color, legend_data$pathway)]
  plot(component_i, vertex.color="white", vertex.label.dist=1, vertex.label.color="black", vertex.label.degree=0, vertex.size=2, main=paste("Component", i))
  legend("bottom", legend=legend_data$pathway, fill=legend_data$color)
}

