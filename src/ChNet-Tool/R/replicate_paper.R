# Load necessary libraries
suppressPackageStartupMessages({
  
  library(ggplot2)
  library(optparse)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(reshape2)
  library(corrplot)
  library(gridExtra)
  library(circlize)
  library(devtools)
  library(igraph)
  library(glmnet)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(plotly)
  library(clusterProfiler)
  # Ensure chNet is installed and loaded
  if (!requireNamespace("chNet", quietly = TRUE)) {
    devtools::install_github("Zhangxf-ccnu/chNet", subdir="pkg")
  }
  library(chNet)
})


read_and_analyze <- function(input_file_1, input_file_2) {
  # Read the first input TSV file into a data frame
  input_1 <- read.table(input_file_1, header = TRUE, sep = "\t")
  
  # Read the second input TSV file into a data frame
  input_2 <- read.table(input_file_2, header = TRUE, sep = "\t")
  
  # Set the gene names as rownames for both data frames
  rownames(input_1) <- input_1$Gene
  rownames(input_2) <- input_2$Gene
  
  # Remove the Gene column from both data frames
  input_1 <- input_1[-1]
  input_2 <- input_2[-1]
  
  # Identify the common genes between the two data frames
  common_genes <- intersect(rownames(input_1), rownames(input_2))
  
  # Subset both data frames to include only the common genes
  input_1 <- input_1[common_genes, , drop = FALSE]
  input_2 <- input_2[common_genes, , drop = FALSE]
  
  # Combine the two data frames column-wise
  combined_data <- cbind(as.matrix(input_1), as.matrix(input_2))
  
  # Create a vector indicating the group for each column in the combined data
  group_vector <- c(rep("CD8 Exhausted T-cells", ncol(input_1)), rep("Macrophages", ncol(input_2)))
  
  # Ensure the rownames of the combined data matrix are the common genes
  rownames(combined_data) <- common_genes
  
  # Transpose the combined data matrix for chNet input
  X_matrix <- t(combined_data)
  
  # Apply the chNet algorithm to the transposed data matrix
  Ch_results <- chNet(X = X_matrix, group = group_vector, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  
  # Return a list containing the prepared data and the analysis results
  return(list(data = list(X = X_matrix, genes = rownames(X_matrix), group = group_vector), results = Ch_results))
}


#' @Visualization

# Differential network graph
# Differential network graph
plot_diff_network <- function(diff_net, diff_gene_names) {
  
  # Calculate gene strengths and determine hub genes using the three-sigma rule
  gene_strengths <- strength(diff_net, mode = "all", weights = NULL)
  mean_strength <- mean(gene_strengths)
  sd_strength <- sd(gene_strengths)
  threshold <- mean_strength + 3 * sd_strength
  hub_genes <- names(which(gene_strengths > threshold))
  
  # Filter out non-connected nodes
  non_isolated_nodes <- V(diff_net)[igraph::degree(diff_net, mode = "all") > 0]
  diff_net <- induced_subgraph(diff_net, non_isolated_nodes)
  
  # If the graph has more than 400 nodes, retain only the top 10% by degree
  if (length(V(diff_net)) > 400) {
    degree_threshold <- quantile(igraph::degree(diff_net, mode = "all"), 0.9)
    high_degree_nodes <- V(diff_net)[igraph::degree(diff_net, mode = "all") >= degree_threshold]
    diff_net <- induced_subgraph(diff_net, high_degree_nodes)
  }
  
  # Calculate vertex sizes based on degrees
  degrees <- igraph::degree(diff_net, mode = "all")
  max_size <- 12  # Increased max size
  min_size <- 4   # Decreased min size
  size_range <- max_size - min_size
  vertex_sizes <- (degrees - min(degrees)) / max(degrees) * size_range + min_size
  vertex_sizes[names(vertex_sizes) %in% hub_genes] <- max_size
  
  # Set vertex shapes and colors
  vertex_shapes <- ifelse(names(V(diff_net)) %in% hub_genes, "square", "circle")
  vertex_colors <- ifelse(names(V(diff_net)) %in% diff_gene_names, "red", "blue")
  
  # Adjust edge weights, emphasizing edges connected to hub genes
  edge_weights <- sapply(E(diff_net), function(e) {
    edge_weight <- E(diff_net)$weight[e]
    if (V(diff_net)[ends(diff_net, e)[1]]$name %in% hub_genes ||
        V(diff_net)[ends(diff_net, e)[2]]$name %in% hub_genes) {
      return(edge_weight * 1.5) 
    } else {
      return(edge_weight)
    }
  })
  edge_weights <- edge_weights / max(edge_weights) * 5
  
  # Plot the differential network graph
  plot(diff_net, layout=layout_with_fr(diff_net),
       vertex.shape=vertex_shapes,
       vertex.color=vertex_colors,
       vertex.size=vertex_sizes,
       edge.width=edge_weights,  
       vertex.label=V(diff_net)$name,
       vertex.label.cex=0.7,  # Slightly smaller label size
       vertex.label.color="black",
       edge.arrow.size=0.5, 
       main = "Differential Network Graph",
       main.cex=1.5,
       main.font=2)
  
  # Add legend to the plot
  legend("topright", legend = c("Diff Gene", "Non-Diff Gene"), pch = c(15, 15), col = c("red", "blue"))
}


# Hub genes boxplot
# Hub genes boxplot
plot_hub_genes_boxplot <- function(expression_data_t, grouping_vector, result, dataset) {
  
  # Calculate gene strengths and determine hub genes using the three-sigma rule
  gene_strengths <- strength(result$Diff.net, weights = E(result$Diff.net)$weight)
  mean_strength <- mean(gene_strengths)
  sd_strength <- sd(gene_strengths)
  threshold <- mean_strength + 3 * sd_strength
  hub_genes <- names(which(gene_strengths > threshold))
  
  # Print identified hub genes
  cat("Hub genes identified:\n", paste(hub_genes, collapse = ", "), "\n")
  
  # Extract expression data for hub genes
  hub_gene_data <- dataset$X[, hub_genes]
  
  # Prepare data frame for plotting
  df_for_plot <- data.frame(
    Expression = as.vector(hub_gene_data),
    Gene = rep(hub_genes, each = nrow(hub_gene_data)),
    Group = rep(dataset$group, times = length(hub_genes))
  )
  
  # Create a boxplot of hub gene expression levels
  ggplot(df_for_plot, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(outlier.size = 1.5, outlier.shape = 16) +  # Customize outliers
    scale_fill_manual(values = c("forestgreen", "tomato")) +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    labs(x = "Gene", y = "Expression Levels", title = "Hub Genes Expression Levels") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
}


# Encapsulated function to handle and save the analysis results as PDFs
replicate_analysis_paper <- function(dataset, condition1, condition2, output_path, pdf_name_suffix) {
  
  # Load the dataset dynamically using the dataset name provided
  eval(parse(text = paste0('data("', dataset, '")')))
  dataset <- get(dataset)
  
  # Perform the chNet analysis on the dataset
  # chNet is a network analysis method that identifies differential networks between groups
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  
  # Transpose the expression data for easier plotting
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  # Extract the names of differentially expressed genes
  diff_gene_names <- names(results$diff.gene[results$diff.gene == 1])
  
  # Plot and save the Differential Network Graph as a PDF
  pdf(file.path(output_path, paste0("Differential_network_graph_", pdf_name_suffix, ".pdf")))
  print(plot_diff_network(results$Diff.net, diff_gene_names))
  dev.off()
  
  # Plot and save the Hub Genes Boxplot as a PDF
  pdf(file.path(output_path, paste0("Hub_genes_boxplot_", pdf_name_suffix, ".pdf")))
  print(plot_hub_genes_boxplot(expression_data_t, grouping_vector, results, dataset))
  dev.off()
}


# Parse command line options
option_list <- list(
  
  make_option(c("--input1"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/Paper_data/GSE13159.AML.tsv", help = "Path to the first input file"),
  make_option(c("--condition1"), type = "character", default = "cancer", help = "Condition for the first input"),
  make_option(c("--condition2"), type = "character", default = "normal", help = "Condition for the second input"),
  make_option(c("--output_path"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/Results/replicate_paper", help = "Path to the output directory"),
  make_option(c("--pdf_name_suffix"), type = "character", default = "GSE13159.AML", help = "Suffix for the PDF files")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Ensure all required arguments are provided
if (is.null(args$input1) || is.null(args$condition1) || is.null(args$condition2) || is.null(args$output_path) || is.null(args$pdf_name_suffix)) {
  print_help(parser)
  stop("Missing arguments. Please specify all required files and parameters.", call. = FALSE)
}

# Execute the analysis and generate PDFs
replicate_analysis_paper(args$input1, args$condition1, args$condition2, args$output_path, args$pdf_name_suffix)
