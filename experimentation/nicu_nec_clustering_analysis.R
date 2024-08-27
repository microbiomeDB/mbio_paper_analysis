## Longitudinal pathway network for NICU-NEC


# Setup
setwd("~/Documents")
library(MicrobiomeDB, quietly = TRUE)
library(igraph, quietly = TRUE)
library(ggplot2)
library(stringr)
source("nicu_nec_analysis/utils.R")

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
  
  shared_pathway_network <- createSharedPathwayNetwork(day_species_collection, day_pathways_collection, method="spearman", corrCoeffThreshold = 0.5, pValueThreshold = 0.05)

  V(shared_pathway_network)$color <- unlist(lapply(V(shared_pathway_network)$name, function(v) {ifelse(v %in% c("Clostridium", "Klebsiella"), "blue", "black")}))
  V(shared_pathway_network)$label.color <- V(shared_pathway_network)$color
  
  igraph::plot.igraph(
    shared_pathway_network,
    arrow.mode=0,
    vertex.label.dist=3,
    vertex.label.degree=0,
    vertex.size=2,
    main=paste("Pathway network,", day),
    edge.width=E(shared_pathway_network)$weight/10,
    layout=layout_nicely(shared_pathway_network)
    # layout=layout_in_circle(shared_pathway_network, order(V(shared_pathway_network)$name))
  )
  return(shared_pathway_network)
})

## Perform clustering and such on the graph list

cluster_results <- data.frame(
  nclusters_med=numeric(),
  nclusters_upperq=numeric(),
  nclusters_lowerq=numeric(),
  max_size_med=numeric(),
  max_size_upperq=numeric(),
  max_size_lowerq=numeric(),
  randomized_nclusters_med=numeric(),
  randomized_nclusters_upperq=numeric(),
  randomized_nclusters_lowerq=numeric(),
  randomized_max_size_med=numeric(),
  randomized_max_size_upperq=numeric(),
  randomized_max_size_lowerq=numeric(),
  stringsAsFactors=FALSE
)

for (g in graph_list) {

  ## Cluster the graph 100 times and take the median of the results
  graph_n_clusters <- numeric()
  graph_max_size <- numeric()

  # Remove isolated nodes and self loops for a more fair randomization test
  g_simple <- simplify(g)
  g_simple_noisolates <- delete.vertices(g_simple , which(degree(g_simple)==0))
  
  for (i in 1:100) {
    clustering <- cluster_infomap(g_simple_noisolates)
    graph_n_clusters <- c(graph_n_clusters, length(clustering))
    graph_max_size <- c(graph_max_size, max(unlist(lapply(seq_along(clustering),function(x){length(clustering[[x]])}))))
    # For components
    # comps <- components(g)
    # graph_n_clusters <- c(graph_n_clusters, comps$no) # number of components
    # graph_max_size <- c(graph_max_size, max(comps$csize))
  }

  # Next we need to randomize the network and perform the clustering again.
  # We'll take 1000 randomized graphs.
  graph_randomized_n_clusters <- numeric()
  graph_randomized_max_size <- numeric()
  for (i in 1:1000) {
    g_random <- rewire(g_simple_noisolates, with=keeping_degseq(loops = TRUE, niter = ecount(g)*10))
    clustering <- cluster_infomap(g_random)
    graph_randomized_n_clusters <- c(graph_randomized_n_clusters, length(clustering))
    graph_randomized_max_size <- c(graph_randomized_max_size, max(unlist(lapply(seq_along(clustering),function(x){length(clustering[[x]])}))))
    # For components
    # comps <- components(g_random)
    # graph_randomized_n_clusters <- c(graph_randomized_n_clusters, comps$no) # number of components
    # graph_randomized_max_size <- c(graph_randomized_max_size, max(comps$csize))
  }

  cluster_results <- rbind(cluster_results, data.frame(
    nclusters_med=median(graph_n_clusters),
    nclusters_upperq=quantile(graph_n_clusters, 0.75),
    nclusters_lowerq=quantile(graph_n_clusters, 0.25),
    max_size_med=median(graph_max_size),
    max_size_upperq=quantile(graph_max_size, 0.75),
    max_size_lowerq=quantile(graph_max_size, 0.25),
    randomized_nclusters_med=median(graph_randomized_n_clusters),
    randomized_nclusters_upperq=quantile(graph_randomized_n_clusters, 0.95),
    randomized_nclusters_lowerq=quantile(graph_randomized_n_clusters, 0.05),
    randomized_max_size_med=median(graph_randomized_max_size),
    randomized_max_size_upperq=quantile(graph_randomized_max_size, 0.95),
    randomized_max_size_lowerq=quantile(graph_randomized_max_size, 0.05),
    stringsAsFactors=FALSE
  ))
}



# Combine diagnosis day with cluster results
gg_data <- cbind(data.frame(diagnosis_day=diagnosis_day), cluster_results)

# Plot
ggplot(gg_data, aes(x=diagnosis_day, y=cluster_results.nclusters)) +
  geom_line()

ggplot(gg_data, aes(x=diagnosis_day, y=cluster_results.max_size)) +
  geom_line()

library(plotly)
fig <- plot_ly(gg_data, x = ~diagnosis_day)
fig <- fig %>% add_trace(y = ~nclusters_med, name = 'NEC samples',mode = 'lines+markers') 
fig <- fig %>% add_trace(y = ~randomized_nclusters_lowerq, name='5% quantile', mode = 'lines')
fig <- fig %>% add_trace(y = ~randomized_nclusters_upperq, name='95% quantile', mode = 'lines', fill="tonexty", color='orange')
fig <- fig %>% add_trace(y = ~randomized_nclusters_med, name = 'Randomized networks, median', mode = 'lines+markers')
fig <- fig %>% layout(title = "Number of clusters across diagnosis day for NEC samples and randomized model",
                      xaxis = list(title = "Diagnosis day",
                                   showline = TRUE,
                                   zeroline = TRUE),
                      yaxis = list(title = "Number of communities",
                                   gridcolor = 'rgb(255,255,255)',
                                   showline = TRUE,
                                   zeroline = TRUE,
                                   rangemode="tozero")
                      )
fig

fig2 <- plot_ly(gg_data, x=~diagnosis_day)
fig2 <- fig2 %>% add_trace(y = ~max_size_med, name="med", mode="lines")
fig2 <- fig2 %>% add_trace(y = ~randomized_max_size_lowerq, name="lower", mode="lines")
fig2 <- fig2 %>% add_trace(y = ~randomized_max_size_upperq, name="upper", mode="lines", fill="tonexty")
fig2 <- fig2 %>% add_trace(y = ~randomized_max_size_med, name="random med", mode="lines")
fig2


## I guess this is what you have to do to save svg?

install.packages('reticulate')
reticulate::install_miniconda()
reticulate::conda_install('r-reticulate', 'python-kaleido')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')

