## Function to return the shared pathway network

createSharedPathwayNetwork <- function(taxaCollection, pathwayCollection, method, corrCoeffThreshold, pValueThreshold) {
  
  # Compute the correlation
  pathway_vs_species <- correlation(taxaCollection, pathwayCollection, method = method)
  
  
  # Filter and extract data.table. IMPORTANT, this filtering here is our definition
  # of "correlated". To create the network, we must turn correlation into a binary
  # thing - either a taxon is correlated with a pathway or it isn't. To make this
  # more clear, when i mean correlated and passing filters, i'll write "correlated".
  corr_stats <- data.table::setDT(pathway_vs_species@statistics@statistics)
  # Assume a negative threshold means <= and a positive one means >=
  if (corrCoeffThreshold < 0) {
    corr_stats_filtered <- corr_stats[corr_stats$correlationCoef <= corrCoeffThreshold & corr_stats$pValue <= pValueThreshold, ]
  } else {
    corr_stats_filtered <- corr_stats[corr_stats$correlationCoef >= corrCoeffThreshold & corr_stats$pValue <= pValueThreshold, ]
  }
  
  
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
  # Note, there's no need to simplify the resulting network because above we
  # ensured that there are no duplicate edges, and also we want to keep all the self-loops
  network <- igraph::graph_from_data_frame(network_df, directed=FALSE)

  
  E(network)$sharedPathway <- network_df$sharedPathway


  # Also create the incidence matrix so we can easily find which pathways are connecting nodes?
  corr_incidence <- reshape(corr_stats_filtered[, c("data1", "data2", "correlationCoef")], direction="wide", idvar="data1", timevar="data2")
  # binarize corr_incidence
  binarized_incidence <- corr_incidence[,-1]
  binarized_incidence[is.na(binarized_incidence)] <- 0
  binarized_incidence[binarized_incidence > 0] <- 1
  binarized_matrix <- data.matrix(binarized_incidence, rownames.force = NA)
  rownames(binarized_matrix) <- corr_incidence$data1
  colnames(binarized_matrix) <- str_replace(colnames(binarized_matrix), "correlationCoef.", "")
  # Reorder rows based on components
  component_membership <- components(network)$membership
  components_node_list <- lapply(seq(1:max(component_membership)), function(i) {
    return(V(network)$name[which(component_membership==i)])
  })
  taxa_order <- unlist(components_node_list)
  taxa_order <- rlist::list.reverse(taxa_order)
  binarized_matrix <- binarized_matrix[taxa_order,,drop=FALSE]
  # Reorder columns based on components
  # connected_components_graphs <- lapply(seq(1:max(component_membership)), function(i) {
  #   g <- induced_subgraph(network, components_node_list[[i]], "create_from_scratch")
  # })
  # pathway_order <- unlist(lapply(connected_components_graphs, function(g) {
  #   g_df <- igraph::as_data_frame(g)
  #   return(unique(g_df$sharedPathway))
  # }))
  # pathway_order <- rlist::list.reverse(pathway_order)
  # binarized_matrix <- binarized_matrix[,pathway_order]

  return(list(network=network, incidence_matrix=binarized_matrix))
}

## Rescale values to be between min_value and max_value
rescale_min_max <- function(x, a=1, b=20, min_value=0, max_value=1) {
  return((b - a) * (x - min_value) / (max_value - min_value) + a)
}


## Plot an igraph and label specific nodes and their neighbors
# This could use some significant reworking to make it more clear... but it's functional for now!
plot_igraph_highlight_taxa <- function(the_graph, cool_taxa, cool_taxa_colors, coords, min_vertex_size=1, max_vertex_size=9, min_abundance=0, max_abundance=1, plot_name="") {
  
  # Find all the neighbors of nodes that are in the cool taxa
  klebsiella_vertices <- V(the_graph)[grep("Klebsiella", V(the_graph)$name)]
  clostridium_vertices <- V(the_graph)[grep("Clostridium", V(the_graph)$name)]

  # Find all neighbors of all klebsiella and clostridium nodes
  klebsiella_neighbors <- unique(unlist(lapply(klebsiella_vertices, function(v) {neighbors(the_graph, v)})))
  clostridium_neighbors <- unique(unlist(lapply(clostridium_vertices, function(v) {neighbors(the_graph, v)})))
  # Remove klebsiella nodes from the list

  V(the_graph)$color <- unlist(lapply(V(the_graph)$name, function(v) {
    if(grepl(cool_taxa[1], v)) {
      return(cool_taxa_colors[1])
    } else if(grepl(cool_taxa[2], v)) {
      return(cool_taxa_colors[2])
    } else {
      return("white")
    }
    }))
  V(the_graph)$label.color <- V(the_graph)$color
  
  # Set the stroke color of the vertices
  V(the_graph)$frame.color <- "gray"
  V(the_graph)$frame.color[clostridium_neighbors] <- cool_taxa_colors[1]
  V(the_graph)$frame.color[klebsiella_neighbors] <- cool_taxa_colors[2]
  pre_coords <- as.matrix(coords[match(V(the_graph)$name, coords$node_name), 1:2])
  igraph::plot.igraph(
    the_graph,
    arrow.mode=0,
    vertex.label="",
    vertex.size=rescale_min_max(sqrt(V(the_graph)$med_abundances/pi), a=min_vertex_size, b=max_vertex_size, min_value = min_abundance, max_value=max_abundance), ## vertex.size maps to radius. Rescale for area
    main=plot_name,
    layout=pre_coords,
    rescale=F,
    xlim=c(-0.8,0.8),
    ylim=c(-0.8,0.8)
  )

  ## Apply labels manually
  #Specify x and y coordinates of labels, adjust outward as desired
  x = pre_coords[,1]*1.3
  y = pre_coords[,2]*1.3

  #create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
  angle = ifelse(atan(-(pre_coords[,1]/pre_coords[,2]))*(180/pi) < 0,  90 + atan(-(pre_coords[,1]/pre_coords[,2]))*(180/pi), 270 + atan(-pre_coords[,1]/pre_coords[,2])*(180/pi))

  #Apply the text labels with a loop with angle as srt
  for (i in 1:length(x)) {
    # Label neighbors as well
    if (i %in% klebsiella_neighbors) {
      text(x=x[i], y=y[i], labels=V(the_graph)$name[i], adj=NULL, pos=NULL, cex=.7, col="gray", srt=angle[i], xpd=T)
    }
    if (i %in% clostridium_neighbors) {
      text(x=x[i], y=y[i], labels=V(the_graph)$name[i], adj=NULL, pos=NULL, cex=.7, col="gray", srt=angle[i], xpd=T)
    }
    # Label interesting taxa
    if(any(unlist(lapply(cool_taxa, function(s) {grepl(s, V(the_graph)$name[i])})))) {
      text(x=x[i], y=y[i], labels=V(the_graph)$name[i], adj=NULL, pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
    }
  }
}