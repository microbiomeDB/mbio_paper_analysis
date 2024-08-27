## Longitudinal pathway network for NICU-NEC

install.packages("miscTools")
install.packages("plotly")
install.packages("rlist")
install.packages("ComplexUpset", type = "source")
install.packages("ggstatsplot")

# Setup
setwd("~/Documents")
library(MicrobiomeDB, quietly = TRUE)
library(igraph, quietly = TRUE)
library(ggplot2)
library(stringr)
library(plotly)
library(ggstatsplot)
library(ComplexUpset)
source("nicu_nec_analysis/utils.R")

cool_taxa <- c("Clostridium", "Klebsiella")
cool_taxa_colors <- c("#B39FD6", "#008300")
intersection_color <- "#5D8CE9"
edge_color <- "#A6AFB5"


results <- createNICUNECSubsettedPathwayNetworks()


## Let's plot them all on the same layout for better comparison
# There's a good intersection here. Make into single weighted network to create
# a combined layout. Then use this layout to plot and highlight nodes in 
# both networks
pre_graph <- results$graph_list[[1]]
almost_graph <- results$graph_list[[2]]
post_graph <- results$graph_list[[3]]
control_graph <- results$graph_list[[4]]
pre_df <- igraph::as_long_data_frame(pre_graph)
almost_df <- igraph::as_long_data_frame(almost_graph)
post_df <- igraph::as_long_data_frame(post_graph)
control_df <- igraph::as_long_data_frame(control_graph)
combined_df <- rbind(pre_df, almost_df, post_df, control_df)
# Simplify combined_df to have one edge per pair of nodes and add a column for number of times seen
combined_df$occurrences <- ave(combined_df$from_name, combined_df$from_name, combined_df$to_name, FUN = length)
# igraph is very particular about how the input df is set up, so let's redo it
combined_df <- combined_df[, c("from_name", "to_name", "occurrences")]

combined_graph <- igraph::graph_from_data_frame(combined_df, directed=FALSE)
E(combined_graph)$weight <- combined_df$occurrences

coords <- as.data.frame(layout_in_circle(combined_graph, order(V(combined_graph)$name)))
coords$node_name <- V(combined_graph)$name


lab.locs <- radian.rescale(x=1:length(V(combined_graph)), direction=-1, start=0)

# Plot settings
# Take the square root of the abundances to get the radius of the nodes
min_abundance_radius <- sqrt(min(V(pre_graph)$med_abundances, V(almost_graph)$med_abundances, V(post_graph)$med_abundances, V(control_graph)$med_abundances)/pi)
max_abundance_radius <- sqrt(max(V(pre_graph)$med_abundances, V(almost_graph)$med_abundances, V(post_graph)$med_abundances, V(control_graph)$med_abundances)/pi)
min_abundance <- min(V(pre_graph)$med_abundances, V(almost_graph)$med_abundances, V(post_graph)$med_abundances, V(control_graph)$med_abundances)
max_abundance <- max(V(pre_graph)$med_abundances, V(almost_graph)$med_abundances, V(post_graph)$med_abundances, V(control_graph)$med_abundances)
min_vertex_size <- 3
max_vertex_size <- 9

# Simplify all the plots to remove self-loops
# Plot pre-diagnosis
plot_igraph_highlight_taxa(simplify(pre_graph), cool_taxa, cool_taxa_colors, intersection_color, coords, min_vertex_size, max_vertex_size, min_abundance_radius, max_abundance_radius, "Pre-diagnosis (-inf, 2)", edge_color)

# Plot almost-diagnosis
plot_igraph_highlight_taxa(simplify(almost_graph), cool_taxa, cool_taxa_colors, intersection_color, coords, min_vertex_size, max_vertex_size, min_abundance_radius, max_abundance_radius, "Just before diagnosis [-2, 0)", edge_color)

# Plot post-diagnosis
plot_igraph_highlight_taxa(simplify(post_graph), cool_taxa, cool_taxa_colors, intersection_color, coords, min_vertex_size, max_vertex_size, min_abundance_radius, max_abundance_radius, "Post-diagnosis [0, inf)", edge_color)

# Plot control
plot_igraph_highlight_taxa(simplify(control_graph), cool_taxa, cool_taxa_colors, intersection_color, coords, min_vertex_size, max_vertex_size, min_abundance_radius, max_abundance_radius, "Control", edge_color)




## Now let's make a legend
# We'll have to make a fake network to create the legend
# Create a network with four nodes and no edges
legend_graph <- igraph::make_empty_graph(3)

# Need to assign appropriate mean abundance values to the legend_graph nodes.

legend_graph <- igraph::set_vertex_attr(legend_graph, "size", value=c(min_vertex_size, min_vertex_size + (max_vertex_size - min_vertex_size)/2, max_vertex_size))

# Now assign locations
# It doesn't really matter where the locations are, since i can move them in illustrator
legend_coords <- matrix(c(0, 0, 0, 0.25, 0, 0.5), ncol=2, byrow=TRUE)

# Plot legend
igraph::plot.igraph(
  legend_graph,
  arrow.mode=0,
  vertex.label=c(min_abundance, (max_abundance - min_abundance)/2, max_abundance),
  vertex.label.dist=5,
  vertex.color="white",
  vertex.label.degree = 0,
  vertex.size=V(legend_graph)$size,
  main="Legend",
  layout=legend_coords,
  rescale=F,
  xlim=c(-0.8,0.8),
  ylim=c(-0.8,0.8)
)




####### STARTING HERE SHOULD BE NEW SCRIPT
## Let's look at the klebsiella pnemoniae neighborhood, but only in the almost_graph
# Find all neighbors of klebsiella pnemoniae
kleb_pnemoniae_ego_graph <- igraph::make_ego_graph(almost_graph, 1, which(V(almost_graph)$name == "Klebsiella pneumoniae"))

# Convert to data frame
kleb_pnemoniae_ego_df <- igraph::as_long_data_frame(kleb_pnemoniae_ego_graph[[1]])

# Split sharedPathway into a vector
kleb_pnemoniae_ego_df$sharedPathway <- strsplit(kleb_pnemoniae_ego_df$sharedPathway, ",, ")

# Plot
igraph::plot.igraph(
  kleb_pnemoniae_ego_graph[[1]],
  arrow.mode=0,
  vertex.label=V(kleb_pnemoniae_ego_graph[[1]])$name,
  vertex.size=2,
  main="Klebsiella pnemoniae neighborhood",
  layout=layout_in_circle(kleb_pnemoniae_ego_graph[[1]], order(V(kleb_pnemoniae_ego_graph[[1]])$name))
)

kp_neighbors <- incident(kleb_pnemoniae_ego_graph[[1]], which(V(kleb_pnemoniae_ego_graph[[1]])$name == "Klebsiella pneumoniae"))



## Box plot of the ages for each day
age_df$term <- factor(age_df$term, levels=c("pre", "almost", "post", "control"))
  
ggbetweenstats(age_df, x=term, y=age)


## Create an upset plot of shared pathways for the "almost" day
almost_shared_pathways_incidence <- t(incidence_matrix_list[[2]])
# Convert binary incidence to boolean
almost_shared_pathways_incidence <- almost_shared_pathways_incidence == 1
kleb_species <- c("Klebsiella pneumoniae", "Klebsiella variicola", "Klebsiella quasipneumoniae")
almost_shared_pathways_incidence <- almost_shared_pathways_incidence[, kleb_species]
names(almost_shared_pathways_incidence) <- colnames(almost_shared_pathways_incidence)
incidence_df <- as.data.frame(almost_shared_pathways_incidence)
# remove rows that are all 0
almost_shared_pathways_incidence <- almost_shared_pathways_incidence[rowSums(almost_shared_pathways_incidence) > 0, ]
upset(
  incidence_df,
  kleb_species,
  keep_empty_groups=TRUE,
  stripes='white'
) +
  ggtitle('Overlap of correlated pathways')




# Sanity check
node <- which(V(almost_graph)$name == kleb_species[3])
ids <- incident(almost_graph,node)
print(E(almost_graph)[ids])
quasi_ps <- unique(unlist(strsplit(E(almost_graph)$sharedPathway[ids], ",, ")))
vari_ps <- unique(unlist(strsplit(E(almost_graph)$sharedPathway[ids], ",, ")))
pne_ps <- unique(unlist(strsplit(E(almost_graph)$sharedPathway[ids], ",, ")))


## Checking other shared pathways of interest
# First, find pathways shared in clostridia in the 'almost' graph
clost_nodes <- V(almost_graph)[grep("Clostridium", V(almost_graph)$name)]$name
# Find pathways for all pairs of clost_nodes
clost_edges <- lapply(seq_len(nrow(almost_df)), function(r, df) {
  df_row <- df[r, ]
  if (df_row$from_name %in% clost_nodes & df_row$to_name %in% clost_nodes) {
    print(paste(df_row$from_name,'--', df_row$to_name, ':', df_row$sharedPathway))
    return(unlist(strsplit(df_row$sharedPathway, ",, ")))
  } else {
    return(NULL)
  }
}, almost_df)

# Found no shared pathways between clostridia nodes

# Next find pathways shared between klebsiella and bifidobacteria, and clostridia and bifidobacteria
kleb_nodes <- V(almost_graph)[grep("Klebsiella", V(almost_graph)$name)]$name
bifido_nodes <- V(almost_graph)[grep("Bifidobacterium", V(almost_graph)$name)]$name

# Klebsiella v bifido edges.
kleb_bifido_edges <- lapply(seq_len(nrow(almost_df)), function(r, df) {
  df_row <- df[r, ]
  if ((df_row$from_name %in% kleb_nodes & df_row$to_name %in% bifido_nodes) | (df_row$from_name %in% bifido_nodes & df_row$to_name %in% kleb_nodes)) {
    print(paste(df_row$from_name,'--', df_row$to_name, ':', df_row$sharedPathway))
    return(unlist(strsplit(df_row$sharedPathway, ",, ")))
  } else {
    return(NULL)
  }
}, almost_df)

# Clostridia v bifido edges
clost_bifido_edges <- lapply(seq_len(nrow(almost_df)), function(r, df) {
  df_row <- df[r, ]
  if ((df_row$from_name %in% clost_nodes & df_row$to_name %in% bifido_nodes) | (df_row$from_name %in% bifido_nodes & df_row$to_name %in% clost_nodes)) {
    print(paste(df_row$from_name,'--', df_row$to_name, ':', df_row$sharedPathway))
    return(unlist(strsplit(df_row$sharedPathway, ",, ")))
  } else {
    return(NULL)
  }
}, almost_df)

