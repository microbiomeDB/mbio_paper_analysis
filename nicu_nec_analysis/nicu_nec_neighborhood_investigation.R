## Investigating neighborhoods of the NICU-NEC shared pathway networks

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

# Grab the network just before NEC diagnosis
almost_graph <- results$graph_list[[2]]

# Let's look at the klebsiella pnemoniae neighborhood
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
