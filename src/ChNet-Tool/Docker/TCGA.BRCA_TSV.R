#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
  library(devtools)
  library(igraph)
  library(glmnet)
  library(chNet)
})

# Helper function to determine if a gene is differentially expressed
is_differentially_expressed <- function(expression_vector_group1, expression_vector_group2) {
  # Applying t-test to see if expression levels are significantly different
  test_result <- t.test(expression_vector_group1, expression_vector_group2)
  return(test_result$p.value < 0.05)  # Returns TRUE if p-value is less than 0.05
}

perform_analysis <- function(dataset_name) {
  data(list = dataset_name)
  dataset <- get(dataset_name)
  
  # Assuming the dataset has expression data in a matrix form with rows as genes and columns as samples
  group_vector <- c(rep("Basal", 95), rep("LumA", 231))  # Modify according to your actual data
  
  result = chNet(X = dataset$X, group = group_vector, subsampling = FALSE, R= 20,
                 lambar = 2.825, parallel = FALSE, nCpus = 4)
  
  # Determine differential expression for each gene
  diff_expressed <- apply(dataset$X, 2, function(gene_expression) {  # 2 for columns as genes
    is_differentially_expressed(
      gene_expression[1:95],  # Basal group samples
      gene_expression[96:326]  # LumA group samples
    )
  })
  
  # Store differential expression results in result object
  result$diff_expressed <- diff_expressed
  return(result)
}

# Updated function to save differential network with "weight condition"
save_diff_network <- function(result, output_path) {
  if (!inherits(result$Diff.net, "igraph")) {
    stop("Diff.net is not an igraph object")
  }
  
  diff_edge_edges <- as.data.frame(as_edgelist(result$Diff.net, names = TRUE))
  diff_edge_edges$weight <- E(result$Diff.net)$weight
  
  # Add a column to indicate differential condition
  diff_edge_edges$condition <- ifelse(result$diff_expressed[diff_edge_edges$V1] | result$diff_expressed[diff_edge_edges$V2],
                                      "differential", "non-differential")
  
  output_file_path <- file.path(output_path, "TCGA.BRCA_network.tsv")
  write_tsv(diff_edge_edges, output_file_path)
  
  cat("Results saved to:", output_file_path, "\n")
}

# Define command line options
option_list <- list(
  make_option(c("--dataset"), type = "character", default = "TCGA.BRCA", help = "Name of the dataset to use"),
  make_option(c("--output_path"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/", help = "Path to the output directory")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Ensure all required arguments are provided
if (is.null(args$dataset) || is.null(args$output_path)) {
  print_help(parser)
  stop("Missing arguments. Please specify the dataset and output path.", call. = FALSE)
}

# Execute the analysis and save the results
results <- perform_analysis(args$dataset)
save_diff_network(results, args$output_path)
