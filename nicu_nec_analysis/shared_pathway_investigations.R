## Specific pathway investigations

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


# Generate shared pathway networks subsetted by relative diagnosis period
results <- createNICUNECSubsettedPathwayNetworks()

almost_graph <- results$graph_list[[2]]

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
