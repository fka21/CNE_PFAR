### SCRIPT METADATA ###
# Title: Phylogenetic Tree Construction and Divergence Time Estimation
# Author: Ferenc Kagan
# Date: 2024.08.27
# Description: This script constructs a phylogenetic tree of selected fish species, estimates divergence times using 
#              the DateLife package, and saves the resulting tree in Newick format.

### LIBRARIES ###
# Install and load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load all required libraries using pacman
pacman::p_load(rotl, tidyverse, phangorn, datelife)

### SET WORKING DIRECTORY ###
setwd("~/Documents/Projects/CNE_PFARs/")

### DATA PROCESSING ###
# Define a vector of species names
fish_species <- c("Betta splendens", "Macropodus opercularis", "Anabas testudineus", 
                  "Danio rerio", "Oryzias latipes", "Oreochromis niloticus", 
                  "Astyanax mexicanus", "Nothobranchius furzeri", "Gasterosteus aculeatus", 
                  "Takifugu rubripes")

# Resolve species names using the Open Tree of Life's Taxonomic Name Resolution Service (TNRS)
resolved_names <- rotl::tnrs_match_names(names = fish_species, context_name = "All life")

### PHYLOGENETIC TREE CONSTRUCTION ###
# Generate a phylogenetic tree using the resolved OTT IDs
phylo_tree <- rotl::tol_induced_subtree(ott_ids = resolved_names$ott_id)

# Plot the initial tree
plot(phylo_tree)

# Clean up node and tip labels
phylo_tree$node.label <- ""
phylo_tree$tip.label <- str_remove_all(phylo_tree$tip.label, "_ott.*")

### DIVERGENCE TIME ESTIMATION ###
# Estimate divergence times using the DateLife package
datelife_results <- datelife::get_datelife_result(phylo_tree)

# Summarize the divergence times into a phylogenetic tree format
dated_tree <- datelife::summarize_datelife_result(datelife_results, summary_format = "phylo_biggest")
dated_tree

# Rename the tip labels for better readability
dated_tree$tip.label <- c("drer", "amex", "onil", "olat", "nfur", "bspl", "mope", "ates", "trub", "gacu")

# Plot the dated phylogenetic tree
plot(dated_tree)

### OUTPUT ###
# Save the final phylogenetic tree in Newick format
write.tree(dated_tree, "output/sp_tree.nwk")
