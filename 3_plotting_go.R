# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("topGO", quietly = TRUE))
  BiocManager::install("topGO")
setwd("phd/campylobacter_project/gwas/250103_final_GWAS/4_go_enrichment/go_analysis/")
library(topGO)
library(dplyr)

#-------------------------------------------------------------------------------
# Currently this uses the standard 

# Define input files
background_file <- "consolidated_pangenome_interpro.tsv"  # Replace with your background file
interest_file <- "consolidated_gene_hits_interscanpro.tsv"  # Replace with your genes of interest file, Just need gene names in column "Gene" in this file

# Load background gene-to-GO mapping
bg_data <- read.delim(background_file, header = TRUE, sep = "\t")
bg_geneID2GO <- strsplit(as.character(bg_data$GO_Terms), split = "\\|")
names(bg_geneID2GO) <- bg_data$Gene
bg_geneID2GO_list <- lapply(bg_geneID2GO, function(x) unique(unlist(x)))

# Load genes of interest
interest_data <- read.delim(interest_file, header = TRUE, sep = "\t")
genes_of_interest <- interest_data$Gene

# Create a gene list for topGO
geneList <- factor(as.integer(names(bg_geneID2GO_list) %in% genes_of_interest))
names(geneList) <- names(bg_geneID2GO_list)

# Function to run topGO for a specific ontology
run_topGO <- function(ontology, geneList, geneID2GO_list) {
  GOdata <- new("topGOdata",
                description = paste("GO analysis for", ontology),
                ontology = ontology,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO_list)
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  topNodes <- GenTable(GOdata, classic = resultFisher, topNodes = 10) # change this based on the type of analysis desired, resultElim or resultFisher
  return(topNodes)
}

# Run and display results for each ontology
cat("\nBiological Process (BP) results:\n")
bp_results <- run_topGO("BP", geneList, bg_geneID2GO_list)
print(bp_results)

cat("\nCellular Component (CC) results:\n")
cc_results <- run_topGO("CC", geneList, bg_geneID2GO_list)
print(cc_results)

cat("\nMolecular Function (MF) results:\n")
mf_results <- run_topGO("MF", geneList, bg_geneID2GO_list)
print(mf_results)

library(ggplot2)

plot_go_terms <- function(results, title) {
  results$log_p_value <- -log10(as.numeric(results$classic))
  
  ggplot(results, aes(x = reorder(Term, -log_p_value), y = log_p_value)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = title, x = "GO Term", y = "-log10(p-value)") +
    theme_minimal()
}

# Generate plots for BP, CC, and MF results
plot_go_terms(bp_results, "Top GO Terms in Biological Process (BP)")
plot_go_terms(cc_results, "Top GO Terms in Cellular Component (CC)")
plot_go_terms(mf_results, "Top GO Terms in Molecular Function (MF)")


# Combine BP, CC, and MF results into a single data frame
combined_results <- rbind(
  cbind(bp_results, Ontology = "BP"),
  cbind(cc_results, Ontology = "CC"),
  cbind(mf_results, Ontology = "MF")
)

# Create log-transformed p-values for plotting
combined_results$log_p_value <- -log10(as.numeric(combined_results$classic))

# Filter to keep only the top 7 GO terms per Ontology based on p-value
top_go_terms <- combined_results %>%
  group_by(Ontology) %>%
  top_n(5, wt = log_p_value) %>%
  ungroup()

# Reorder bars within each GO category by p-value
top_go_terms$Term <- factor(top_go_terms$Term, 
                            levels = top_go_terms$Term[order(top_go_terms$log_p_value, decreasing = TRUE)])
library(ggpubr)
# Plot the data
ggplot(top_go_terms, aes(x = Term, y = log_p_value, fill = Ontology)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~Ontology, scales = "free_y") +  # Separate by GO category, each with its own y-axis scale
  labs(title = "Top 5 GO Terms in Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)",
       x = "GO Term", y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide legend, since color is already used for facets
