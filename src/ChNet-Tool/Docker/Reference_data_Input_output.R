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

# Define the main function for performing analysis
perform_analysis <- function(input_file_1, input_file_2) {
  # Load the datasets using correct paths
  input_1 <- read.table(input_file_1, header = TRUE, sep = "\t")
  input_2 <- read.table(input_file_2, header = TRUE, sep = "\t")
  
  # Convert tibbles to data frames and set row names
  rownames(input_1) <- input_1$Gene
  input_1 <- input_1[,-1]
  rownames(input_2) <- input_2$Gene
  input_2 <- input_2[,-1]
  
  # Ensure that both datasets have the same genes in the same order
  all_genes <- intersect(rownames(input_1), rownames(input_2))
  
  input_1 <- input_1[all_genes, , drop=FALSE]
  input_2 <- input_2[all_genes, , drop=FALSE]
  
  # Combine the two datasets
  combined_data <- cbind(as.matrix(input_1), as.matrix(input_2))
  group_vector <- c(rep("CD8_exhausted", ncol(input_1)), rep("Macrophages", ncol(input_2)))
  
  rownames(combined_data) <- all_genes
  X_matrix <- t(combined_data)  # Transpose so that genes are columns
  
  # Run analysis
  chNet_results <- chNet(X = X_matrix, group = group_vector, subsampling = FALSE, R = 20, lambar = 2.8, parallel = FALSE, nCpus = 4)
  
  # Check if Diff.net is an igraph object
  if (!inherits(chNet_results$Diff.net, "igraph")) {
    stop("Expected Diff.net to be an igraph object, check chNet function output")
  }
  
  return(list(data = list(X = X_matrix, genes = rownames(X_matrix), group = group_vector), results = chNet_results))
}

# Function to save differential network results
save_diff_network <- function(result, output_dir, file_name) {
  if (!inherits(result$results$Diff.net, "igraph")) {
    stop("Diff.net is not an igraph object")
  }
  
  diff_edge_edges <- as.data.frame(as_edgelist(result$results$Diff.net, names = TRUE))
  diff_edge_edges$weight <- E(result$results$Diff.net)$weight
  
  diff_status <- result$results$diff.gene
  
  label_edge_condition <- function(edge_from, edge_to, diff_status) {
    condition <- ifelse(diff_status[edge_from] == 1 && diff_status[edge_to] == 1, "differential", "non-differential")
    return(condition)
  }
  
  diff_edge_edges$condition <- mapply(label_edge_condition,
                                      diff_edge_edges$V1,
                                      diff_edge_edges$V2,
                                      MoreArgs = list(diff_status = diff_status))
  
  diff_edge_edges <- diff_edge_edges[, c("V1", "V2", "condition", "weight")]
  names(diff_edge_edges) <- c("regulator", "target", "condition", "weight")
  
  # Save the formatted results as a TSV file
  output_file_path <- file.path(output_dir, paste0("network_ref_", file_name, ".tsv"))
  write.table(diff_edge_edges, output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Print path to output file
  cat("Results saved to:", output_file_path, "\n")
}

# Handling the reference datasets
chNet_references_data <- function(input1, input2, output_dir, file_name_suffix) {
  ref_output <- perform_analysis(input1, input2)
  save_diff_network(ref_output, output_dir, file_name_suffix)
}

# Define command line options
option_list <- list(
  make_option(c("--input1"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/Ref_dataset/500/out_CD8_exhausted.tsv", metavar = "FILE"),
  make_option(c("--input2"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/Ref_dataset/500/out_Macrophages.tsv", metavar = "FILE"),
  make_option(c("--output_path"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/", metavar = "DIR"),
  make_option(c("--file_name_suffix"), type = "character", default = "500", help = "Suffix for the files")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Ensure all required arguments are provided
if (args$input1 == "" || args$input2 == "" || args$output_path == "" || args$file_name_suffix == "") {
  print_help(parser)
  stop("Missing arguments. Please specify all required files and output path.", call. = FALSE)
}

# Run the reference data analysis
chNet_references_data(args$input1, args$input2, args$output_path, args$file_name_suffix)
