# Load necessary libraries
suppressPackageStartupMessages({

library(chNet)
library(ggplot2)
library(optparse)
library(dplyr)
library(readr)
library(RColorBrewer)
library(reshape2)
library(corrplot)
library(gridExtra)
library(circlize)
library(igraph)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(plotly)
library(clusterProfiler)
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


#' @DownStream Visualization
#Clustering Plot
network_clustering <- function(diff_net) {
  # Calculate gene strengths and identify hub genes using the three-sigma rule
  gene_strengths <- strength(diff_net, mode = "all", weights = NULL)
  mean_strength <- mean(gene_strengths)
  sd_strength <- sd(gene_strengths)
  threshold <- mean_strength + 3 * sd_strength
  hub_genes <- names(which(gene_strengths > threshold))
  
  # Remove isolated nodes
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
  max_size <- 13
  min_size <- 5
  size_range <- max_size - min_size
  vertex_sizes <- (degrees - min(degrees)) / max(degrees) * size_range + min_size
  vertex_sizes[names(vertex_sizes) %in% hub_genes] <- max_size
  
  # Set vertex shapes
  vertex_shapes <- ifelse(names(V(diff_net)) %in% hub_genes, "square", "circle")
  
  # Determine community structure
  communities <- cluster_louvain(diff_net)
  vertex_colors <- viridis::viridis(max(communities$membership))[communities$membership]
  
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
       vertex.label.cex=0.8,
       vertex.label.color="red",
       edge.color="gray50",
       edge.arrow.size=0.5, 
       main = "Clustering of Differential Network Graph",
       main.cex=1.5,
       main.font=2,
       sub = "Gene interaction network with hub genes highlighted",
       sub.cex = 1.2,
       sub.font = 3)
}



## Boxplot for the most significant genes
hub_gene_boxplot <- function(dataset, hub_gene_data, top_hub_genes){
  # Prepare data frame for plotting
  df_for_plot <- data.frame(
    Expression = as.vector(hub_gene_data),
    Gene = rep(top_hub_genes, each = nrow(hub_gene_data)),
    Group = rep(dataset$group, times = length(top_hub_genes))
  )
  
  # Create a boxplot of gene expression levels
  ggplot(df_for_plot, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2, notch = TRUE) +
    scale_fill_manual(values = c("forestgreen", "tomato")) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "top"
    ) +
    labs(x = "Gene", y = "Expression Levels", title = "Hub Gene Expression Boxplot")
}


#Correlation for the Most Significant Genes
hub_gene_correlation <- function(dataset, hub_gene_data, top_hub_genes){
  # Calculate the correlation matrix
  cor_matrix <- cor(hub_gene_data)
  
  # Create a correlation plot
  corrplot(cor_matrix, method = "circle", type = "upper",
           order = "hclust",
           tl.col = "darkblue",
           tl.srt = 45,
           tl.cex = 0.8,
           addCoef.col = "darkred",
           number.cex = 0.7,
           main = "Hub Gene Expression Correlation Plot",
           mar = c(4, 4, 2, 2),
           cl.pos = "b")
}



# Gene Ontology Enrichment - Molecular Functions of hub genes  
Gene_Ontology_Enrichment_MF <- function(dataset, numHubGenes) {
  # Run chNet to get differential network results
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  diff_net <- results$Diff.net
  
  # Calculate gene strengths and select top hub genes
  gene_strengths <- strength(diff_net, weights = E(diff_net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:numHubGenes])
  
  # Map gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = top_hub_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  if (length(entrez_ids) == 0) {
    cat("No Entrez IDs could be retrieved for the top hub genes.\n")
    return(NULL)
  }
  
  # Perform GO enrichment analysis for Molecular Functions
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = 'ENTREZID',
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.1,
                  qvalueCutoff = 0.1,
                  readable = TRUE)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    df <- as.data.frame(ego)
    df <- df[order(df$pvalue), ]
    top_df <- head(df, numHubGenes)
    if (nrow(top_df) > 0) {
      # Combine descriptions with gene names for the x-axis
      top_df$label <- paste(top_df$Description, "\n(Genes: ", top_df$geneID, ")", sep = "")
      
      # Create a bar plot for the top GO terms
      p <- ggplot(top_df, aes(x = reorder(label, pvalue), y = Count)) +
        geom_col(fill = "yellow") +
        geom_text(aes(label = Count), vjust = -0.3, color = "black", size = 4) +
        labs(x = "Molecular Functions and Genes", y = "Number of Connected Genes", title = "Top Hub Genes GO Enrichment Analysis - Molecular Functions") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 75, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "none"
        )
      
      return(p)
    } else {
      cat("No significant GO terms identified for the MF ontology.\n")
      return(NULL)
    }
  } else {
    cat("No significant GO terms identified for the MF ontology.\n")
    return(NULL)
  }
}


# Gene Ontology Enrichment - Biological Process of hub genes  
Gene_Ontology_Enrichment_BP <- function(dataset, numHubGenes) {
  # Run chNet to get differential network results
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  diff_net <- results$Diff.net
  
  # Calculate gene strengths and select top hub genes
  gene_strengths <- strength(diff_net, weights = E(diff_net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:numHubGenes])
  
  # Map gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = top_hub_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  if (length(entrez_ids) == 0) {
    cat("No Entrez IDs could be retrieved for the top hub genes.\n")
    return(NULL)
  }
  
  # Perform GO enrichment analysis for Biological Processes
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = 'ENTREZID',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.1,
                  qvalueCutoff = 0.1,
                  readable = TRUE)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    df <- as.data.frame(ego)
    df <- df[order(df$pvalue), ]
    top_df <- head(df, numHubGenes)
    if (nrow(top_df) > 0) {
      # Combine descriptions with gene names for the x-axis
      top_df$label <- paste(top_df$Description, "\n(Genes: ", top_df$geneID, ")", sep = "")
      
      # Create a bar plot for the top GO terms
      p <- ggplot(top_df, aes(x = reorder(label, pvalue), y = Count)) +
        geom_col(fill = "yellow") +
        geom_text(aes(label = Count), vjust = -0.3, color = "black", size = 4) +
        labs(x = "Biological Process and Genes", y = "Number of Connected Genes", title = "Top Hub Genes GO Enrichment Analysis - Biological Process") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 75, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "none"
        )
      
      return(p)
    } else {
      cat("No significant GO terms identified for the BP ontology.\n")
      return(NULL)
    }
  } else {
    cat("No significant GO terms identified for the BP ontology.\n")
    return(NULL)
  }
}



#' @references_data
# Function to handle and analyze reference datasets
downstream_analysis_ref <- function(input1, input2, condition1, condition2, output_path, pdf_name_suffix) {
  # Read and analyze the input files
  ref_output <- read_and_analyze(input1, input2)
  dataset <- ref_output$data
  results <- ref_output$results
  
  # Transpose expression data and get grouping vector
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  # Extract the differential network
  diff_net <- results$Diff.net
  
  # Calculate gene strengths and identify the top 10 hub genes
  gene_strengths <- strength(results$Diff.net, weights = E(results$Diff.net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:10])
  hub_gene_data <- dataset$X[, top_hub_genes]
  
  # Generate and save the network clustering plot
  pdf(file.path(output_path, paste0("network_clustering_", pdf_name_suffix, ".pdf")))
  result <- network_clustering(diff_net)
  if (!is.null(result$plot)) { # Avoid unnecessary print
    print(result$plot)
  }
  dev.off()
  
  # Generate and save the hub gene expression boxplot
  pdf(file.path(output_path, paste0("hub_gene_boxplot_", pdf_name_suffix, ".pdf")))
  print(hub_gene_boxplot(dataset, hub_gene_data, top_hub_genes))
  dev.off()
  
  # Generate and save the hub gene correlation plot
  pdf(file.path(output_path, paste0("hub_gene_correlation_", pdf_name_suffix, ".pdf")))
  hub_gene_correlation(dataset, hub_gene_data, top_hub_genes)
  dev.off()
  
  # Generate and save the GO enrichment plot for Biological Process
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_BP_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_BP(dataset, numHubGenes = 10))
  dev.off()
  
  # Generate and save the GO enrichment plot for Molecular Function
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_MF_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_MF(dataset, numHubGenes = 10))
  dev.off()
}


#' @Paper_dataset
# Function to handle and analyze datasets from the paper
downstream_analysis_paper <- function(dataset, condition1, condition2, output_path, pdf_name_suffix) {
  # Load the dataset
  eval(parse(text = paste0('data("', dataset, '")')))
  dataset <- get(dataset)
  
  # Apply the chNet algorithm
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R= 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  
  # Transpose expression data and get grouping vector
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  # Extract the differential network
  diff_net <- results$Diff.net
  
  # Calculate gene strengths and identify the top 10 hub genes
  gene_strengths <- strength(results$Diff.net, weights = E(results$Diff.net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:10])
  hub_gene_data <- dataset$X[, top_hub_genes]
  
  # Generate and save the GO enrichment plot for Biological Process
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_BP_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_BP(dataset, numHubGenes = 10))
  dev.off()
  
  # Generate and save the GO enrichment plot for Molecular Function
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_MF_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_MF(dataset, numHubGenes = 10))
  dev.off()
  
  # Generate and save the network clustering plot
  pdf(file.path(output_path, paste0("network_clustering_", pdf_name_suffix, ".pdf")))
  result <- network_clustering(diff_net)
  if (!is.null(result$plot)) { # Avoid unnecessary print
    print(result$plot)
  }
  dev.off()
  
  # Generate and save the hub gene expression boxplot
  pdf(file.path(output_path, paste0("hub_gene_boxplot_", pdf_name_suffix, ".pdf")))
  print(hub_gene_boxplot(dataset, hub_gene_data, top_hub_genes))
  dev.off()
  
  # Generate and save the hub gene correlation plot
  pdf(file.path(output_path, paste0("hub_gene_correlation_", pdf_name_suffix, ".pdf")))
  hub_gene_correlation(dataset, hub_gene_data, top_hub_genes)
  dev.off()
}



# Parse command line options
option_list <- list(
  make_option(c("--operation_type"), type = "character", default = "ref", help = "paper or ref"),
  make_option(c("--input1"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/ref_dataset/500/out_CD8_exhausted.tsv", help = "Path to the first input file"),
  make_option(c("--input2"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/ref_dataset/500/out_Macrophages.tsv",help = "Path to the second input file"),
  make_option(c("--condition1"), type = "character",default = "CD8 Exhausted T-cells" , help = "Condition for the first input"),
  make_option(c("--condition2"), type = "character",default = "Macrophages" , help = "Condition for the second input"),
  make_option(c("--output_path"), type = "character", default = "C:/Users/Piyal/Desktop/Bio/Results/downstream/500", help = "Path to the output directory"),
  make_option(c("--pdf_name_suffix"), type = "character",default = "500", help = "Suffix for the PDF files")
)

  
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Ensure all required arguments are provided
if (is.null(args$operation_type) || is.null(args$input1) || is.null(args$condition1) || is.null(args$condition2) || is.null(args$output_path) || is.null(args$pdf_name_suffix) || (args$operation_type == "ref" && is.null(args$input2))) {
  print_help(parser)
  stop("Missing arguments. Please specify all required files and parameters.", call. = FALSE)
}

# Call the appropriate function based on the operation type
if (!is.null(args$operation_type) && args$operation_type == "ref") {
  downstream_analysis_ref(args$input1, args$input2, args$condition1, args$condition2, args$output_path, args$pdf_name_suffix)
} else if (!is.null(args$operation_type) && args$operation_type == "paper") {
  downstream_analysis_paper(args$input1, args$condition1, args$condition2, args$output_path, args$pdf_name_suffix)
} else {
  stop("Invalid operation type. Please specify 'paper' or 'ref'.")
}
