## Sanity chekcks

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


# Sanity check
node <- which(V(almost_graph)$name == kleb_species[3])
ids <- incident(almost_graph,node)
print(E(almost_graph)[ids])
quasi_ps <- unique(unlist(strsplit(E(almost_graph)$sharedPathway[ids], ",, ")))
vari_ps <- unique(unlist(strsplit(E(almost_graph)$sharedPathway[ids], ",, ")))
pne_ps <- unique(unlist(strsplit(E(almost_graph)$sharedPathway[ids], ",, ")))

