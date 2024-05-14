#!/usr/bin/env Rscript

organism = "org.Dm.eg.db"
# Load necessary libraries
suppressPackageStartupMessages({
library(organism, character.only = TRUE)
library(readr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(optparse)
library(pbapply)
})

# Set up command-line options
option_list <- list(
  make_option(c("-f", "--network_file"), type = "character", default = NULL, help = "Path to the TSV file containing network data"),
  make_option(c("-o", "--output_path"), type = "character", default = NULL, help = "Output directory to save the enrichment map image"),
  make_option(c("--ont"), type = "character", default = "ALL", help = "Ontology to use for GSEA. Options: 'BP', 'CC', 'MF', 'ALL' (Biological Process, Cellular Component, Molecular Function, or All)"),
  make_option(c("--keyType"), type = "character", default = "ALIAS", help = "Key type for GSEA. Common options: 'ENTREZID', 'SYMBOL', 'ALIAS'"),
  make_option(c("--nPerm"), type = "integer", default = 10000, help = "Number of permutations for GSEA, determines the accuracy of the result"),
  make_option(c("--minGSSize"), type = "integer", default = 3, help = "Minimum gene set size for GSEA"),
  make_option(c("--maxGSSize"), type = "integer", default = 800, help = "Maximum gene set size for GSEA"),
  make_option(c("--pvalueCutoff"), type = "double", default = 0.05, help = "P-value cutoff for GSEA"),
  make_option(c("--verbose"), type = "logical", default = TRUE, help = "Verbose output. Options: TRUE or FALSE"),
  make_option(c("--pAdjustMethod"), type = "character", default = "none", help = "P-value adjustment method for GSEA. Options include 'BH', 'BY', 'holm', 'Bonferroni', 'none'")
)

# Parse options
args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args$network_file) || is.null(args$output_path)) {
  stop("Both --network_file and --output_path must be provided.", call. = FALSE)
}

# Load the network data and provide feedback
message("Loading network data from: ", args$network_file)
network_data <- read_tsv(args$network_file)

plot_gsea <- function(data, weight_column, output_filename, plot_title, args) {
  message("Preparing gene list for GSEA...")
  gene_list <- data %>%
    select(Target, Regulator, !!sym(weight_column)) %>%
    pivot_longer(cols = c(Target, Regulator), names_to = "Gene_role", values_to = "Gene") %>%
    group_by(Gene) %>%
    summarise(Score = mean(!!sym(weight_column), na.rm = TRUE), .groups = 'drop')
  
  message("Performing Gene Set Enrichment Analysis...")
  # Convert the data frame to a named vector
  gene_list_vector <- setNames(gene_list$Score, gene_list$Gene)
  
  # Sort the gene list by z-score in decreasing order
  gene_list_sorted <- sort(gene_list_vector, decreasing = TRUE)
  
  # Apply GSE
  gse <- gseGO(
    geneList = gene_list_sorted,
    ont = args$ont,
    keyType = args$keyType,
    nPerm = args$nPerm,
    minGSSize = args$minGSSize,
    maxGSSize = args$maxGSSize,
    pvalueCutoff = args$pvalueCutoff,
    verbose = args$verbose,
    OrgDb = organism,
    pAdjustMethod = args$pAdjustMethod
  )

  
  if (length(gse) > 0) {
    message("Plotting enrichment map...")
    gse_sim <- pairwise_termsim(gse)
    p <- emapplot(gse_sim, showCategory = 10) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    ggsave(filename = file.path(args$output_path, output_filename), plot = p, width = 10, height = 8)
    message("Enrichment map saved to: ", file.path(args$output_path, output_filename))
  } else {
    message("No significant pathways were enriched based on the input gene list.")
  }
}

plot_condition1_gsea <- function(network_data, args) {
  filtered_data <- network_data %>% filter(Condition == "condition1")
  plot_gsea(filtered_data, "Weight_1", "condition_1_enrichment_map.png", "Condition 1 Enrichment Map", args)
}

plot_condition2_gsea <- function(network_data, args) {
  filtered_data <- network_data %>% filter(Condition == "condition2")
  plot_gsea(filtered_data, "Weight_1", "condition_2_enrichment_map.png", "Condition 2 Enrichment Map", args)
}


# Execute functions with progress indication
message("Starting combine analysis...")
plot_gsea(network_data, "Weight_3", "diff_corr_enrichment_map.png", "Differentially Correlated Enrichment Map of Gene Sets", args)
message("Analysis completed.")

message("Starting condition 1 analysis...")
plot_condition1_gsea(network_data, args)
message("Analysis completed.")

message("Starting condition 2 analysis...")
plot_condition2_gsea(network_data, args)
message("Analysis completed.")