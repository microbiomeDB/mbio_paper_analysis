## Supporting plots for the NICU NEC analysis

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

age_df <- results$age_df

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

