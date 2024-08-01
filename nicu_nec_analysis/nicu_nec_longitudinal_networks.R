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
  "Shotgun metagenomics Species (Relative taxonomic abundance analysis)",
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
diagnosis_day <- c("pre", "almost", "post", "control") # For making ranges
cool_taxa <- c("Clostridium", "Klebsiella")


# Create a list of correlation graphs, one for each age in ages.
graph_list <- lapply(diagnosis_day, function(day) {
  
  # Subset to the appropriate abundances for this age
  if (day=="pre") {
    day_samples <- sampleMetadata$Sample_Id[!is.na(sampleMetadata$days_of_period_nec_diagnosed_days) & sampleMetadata$days_of_period_nec_diagnosed_days < -2]
  }
  else if (day=="almost") {
    day_samples <- sampleMetadata$Sample_Id[!is.na(sampleMetadata$days_of_period_nec_diagnosed_days) & sampleMetadata$days_of_period_nec_diagnosed_days >= -2 & sampleMetadata$days_of_period_nec_diagnosed_days <0]
  } else if (day=="post") {
    day_samples <- sampleMetadata$Sample_Id[!is.na(sampleMetadata$days_of_period_nec_diagnosed_days) & sampleMetadata$days_of_period_nec_diagnosed_days >= 0]
  } else {
    day_samples <- sampleMetadata$Sample_Id[is.na(sampleMetadata$days_of_period_nec_diagnosed_days)]
  }
  
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
  
  V(simplified_network)$color <- unlist(lapply(V(simplified_network)$name, function(v) {ifelse(v %in% cool_taxa, "blue", "black")}))
  V(simplified_network)$label.color <- V(simplified_network)$color
  
  igraph::plot.igraph(
    simplified_network,
    arrow.mode=0,
    vertex.label.dist=3,
    vertex.label.degree=0,
    vertex.size=2,
    main=paste("Pathway network,", day),
    edge.width=E(simplified_network)$weight/10,
    layout=layout_in_circle(simplified_network)
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


## Let's plot them all on the same layout for better comparison
# There's a good intersection here. Make into single weighted network to create
# a combined layout. Then use this layout to plot and highlight nodes in 
# both networks
pre_df <- igraph::as_long_data_frame(graph_list[[1]])
almost_df <- igraph::as_long_data_frame(graph_list[[2]])
post_df <- igraph::as_long_data_frame(graph_list[[3]])
control_df <- igraph::as_long_data_frame(graph_list[[4]])
combined_df <- rbind(pre_df, almost_df, post_df, control_df)
# Simplify combined_df to have one edge per pair of nodes and add a column for number of times seen
combined_df$occurrences <- ave(combined_df$from_name, combined_df$from_name, combined_df$to_name, FUN = length)
# igraph is very particular about how the input df is set up, so let's redo it
combined_df <- combined_df[, c("from_name", "to_name", "occurrences")]

combined_graph <- igraph::graph_from_data_frame(combined_df, directed=FALSE)
E(combined_graph)$weight <- combined_df$occurrences
pre_graph <- graph_list[[1]]
almost_graph <- graph_list[[2]]
post_graph <- graph_list[[3]]
control_graph <- graph_list[[4]]
coords <- as.data.frame(layout_in_circle(combined_graph))
coords$node_name <- V(combined_graph)$name

## (from stack overflow) Get the labels aligned consistently around the edge of the circle
## for any n of nodes.
## This code borrows bits of ggplot2's polar_coord function
## start = offset from 12 o'clock in radians
## direction = 1 for clockwise; -1 for anti-clockwise.

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
lab.locs <- radian.rescale(x=1:length(V(combined_graph)), direction=-1, start=0)

# Plot pre-diagnosis
V(pre_graph)$color <- unlist(lapply(V(pre_graph)$name, function(v) {ifelse(v %in% cool_taxa, "blue", "black")}))
V(pre_graph)$label.color <- V(pre_graph)$color
pre_coords <- as.matrix(coords[match(V(pre_graph)$name, coords$node_name), 1:2])
igraph::plot.igraph(
  pre_graph,
  arrow.mode=0,
  vertex.label="",
  # vertex.label.dist=4,
  # vertex.label.degree=lab.locs[which(V(pre_graph)$name %in% V(combined_graph)$name)],
  vertex.size=5,
  main="pre-diagnosis",
  layout=pre_coords,
  rescale=F,xlim=c(-0.8,0.8),ylim=c(-0.8,0.8)
)

## Apply labels manually
#Specify x and y coordinates of labels, adjust outward as desired
x = pre_coords[,1]*1.3
y = pre_coords[,2]*1.3

#create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(pre_coords[,1]/pre_coords[,2]))*(180/pi) < 0,  90 + atan(-(pre_coords[,1]/pre_coords[,2]))*(180/pi), 270 + atan(-pre_coords[,1]/pre_coords[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(pre_graph)$name[i], adj=NULL, pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
}

## Almost graph
V(almost_graph)$color <- unlist(lapply(V(almost_graph)$name, function(v) {ifelse(v %in% cool_taxa, "blue", "black")}))
V(almost_graph)$label.color <- V(almost_graph)$color
almost_coords <- as.matrix(coords[match(V(almost_graph)$name, coords$node_name), 1:2])
igraph::plot.igraph(
  almost_graph,
  arrow.mode=0,
  vertex.label.dist=1,
  vertex.label.degree=0,
  vertex.size=10,
  main="just-before-diagnosis",
  layout=almost_coords,
  rescale=F,xlim=c(-0.8,0.8),ylim=c(-0.8,0.8)
)
## Apply labels manually
#Specify x and y coordinates of labels, adjust outward as desired
subgraph_coords <- post_coords
x = subgraph_coords[,1]*1.3
y = subgraph_coords[,2]*1.3

#create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(subgraph_coords[,1]/subgraph_coords[,2]))*(180/pi) < 0,  90 + atan(-(subgraph_coords[,1]/subgraph_coords[,2]))*(180/pi), 270 + atan(-subgraph_coords[,1]/subgraph_coords[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(post_graph)$name[i], adj=NULL, pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
}





V(post_graph)$color <- unlist(lapply(V(post_graph)$name, function(v) {ifelse(v %in% cool_taxa, "blue", "black")}))
V(post_graph)$label.color <- V(post_graph)$color
post_coords <- as.matrix(coords[match(V(post_graph)$name, coords$node_name), 1:2])
igraph::plot.igraph(
  post_graph,
  arrow.mode=0,
  vertex.label.dist=1,
  vertex.label.degree=0,
  vertex.size=10,
  main='post-diagnosis',
  layout=post_coords,
  rescale=F,xlim=c(-0.8,0.8),ylim=c(-0.8,0.8)
)


V(control_graph)$color <- unlist(lapply(V(control_graph)$name, function(v) {ifelse(v %in% cool_taxa, "blue", "black")}))
V(control_graph)$label.color <- V(control_graph)$color
control_coords <- as.matrix(coords[match(V(control_graph)$name, coords$node_name), 1:2])
igraph::plot.igraph(
  control_graph,
  arrow.mode=0,
  vertex.label.dist=1,
  vertex.label.degree=0,
  vertex.size=10,
  main='control',
  layout=control_coords,
  rescale=F,xlim=c(-0.8,0.8),ylim=c(-0.8,0.8)
)






