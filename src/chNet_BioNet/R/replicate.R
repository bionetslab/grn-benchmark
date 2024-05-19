library(chNet)
library(ggplot2)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(pheatmap)


#' @preprocessing
# read and merge two tsv files and apply the chNet algorithm
read_and_analyze <- function(file_path1, file_path2) {
  data1 <- read.table(file_path1, header = TRUE, sep = "\t")
  data2 <- read.table(file_path2, header = TRUE, sep = "\t")
  rownames(data1) <- data1$Gene
  rownames(data2) <- data2$Gene
  data1 <- data1[-1] 
  data2 <- data2[-1]
  
  common_genes <- intersect(rownames(data1), rownames(data2))
  
  data1 <- data1[common_genes, , drop = FALSE]
  data2 <- data2[common_genes, , drop = FALSE]
  
  combined_data <- cbind(as.matrix(data1), as.matrix(data2))
  group_vector <- c(rep("CD8 Exhausted T-cells", ncol(data1)), rep("Macrophages", ncol(data2)))
  
  rownames(combined_data) <- common_genes
  X_matrix <- t(combined_data)
  analysis_results <- chNet(X = X_matrix, group = group_vector, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  
  return(list(data = list(X = X_matrix, genes = rownames(X_matrix), group = group_vector), results = analysis_results))
}


#' @Visualization
plot_individual_heatmap <- function(expression_data, group_vector, group_name) {
  color_palette <- colorRampPalette(c("green", "black", "red"))(299)
  
  cond_indices <- which(group_vector == group_name)
  data_group_name <- expression_data[, cond_indices, drop = FALSE]
  data_group_name_t <- t(data_group_name)
  
  avg_expression_samples <- colMeans(data_group_name_t)
  sample_categories <- cut(avg_expression_samples, 
                           breaks = quantile(avg_expression_samples, probs = c(0, 0.33, 0.66, 1)), 
                           include.lowest = TRUE, labels = c("Low", "Medium", "High"))
  
  sample_annotations <- HeatmapAnnotation(df = data.frame(sample_category = sample_categories),
                                          col = list(sample_category = c("Low" = "green", "Medium" = "black", "High" = "red")))
  
  hm <- Heatmap(data_group_name_t,
                name = "Expression Level",
                col = color_palette,
                top_annotation = sample_annotations,
                show_row_names = FALSE,
                show_column_names = FALSE,
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "complete",
                clustering_distance_columns = "euclidean",
                clustering_method_columns = "complete",
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                width = unit(50, "mm"),
                row_title = paste0("Gene Expression Data of ", group_name))
  
  draw(hm, padding = unit(2, "mm"))
}


# Difference of Partial Correlation
plot_partial_cor_diff_matrix <- function(diff_edge_weight) {
  color_palette <- colorRampPalette(c("blue", "white", "red"))(200)
  partial_cor_diff_matrix <- as.matrix(diff_edge_weight)
  diag(partial_cor_diff_matrix) <- 0  # Ensure zero diagonal
  pheatmap(partial_cor_diff_matrix, color = color_palette, clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean", cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE, main="Difference of Partial Correlation under two conditions." )
}


# Differential network graph
plot_diff_network <- function(diff_net, diff_gene_names) {

  # Calculate gene strengths and determine hub genes using the three-sigma rule
  gene_strengths <- strength(diff_net, mode = "all", weights = NULL)
  mean_strength <- mean(gene_strengths)
  sd_strength <- sd(gene_strengths)
  threshold <- mean_strength + 3 * sd_strength
  hub_genes <- names(which(gene_strengths > threshold))
  # filtering the nodes
  non_isolated_nodes <- V(diff_net)[igraph::degree(diff_net, mode = "all") > 0]
  diff_net <- induced_subgraph(diff_net, non_isolated_nodes)

  if (length(V(diff_net)) > 400) {
    degree_threshold <- quantile(igraph::degree(diff_net, mode = "all"), 0.9)
    high_degree_nodes <- V(diff_net)[igraph::degree(diff_net, mode = "all") >= degree_threshold]
    diff_net <- induced_subgraph(diff_net, high_degree_nodes)
  }

  degrees <- igraph::degree(diff_net, mode = "all")
  
  max_size <- 11 
  min_size <- 5  
  size_range <- max_size - min_size
  vertex_sizes <- (degrees - min(degrees)) / max(degrees) * size_range + min_size
  vertex_sizes[names(vertex_sizes) %in% hub_genes] <- max_size
  vertex_shapes <- ifelse(names(V(diff_net)) %in% hub_genes, "square", "circle")
  vertex_colors <- ifelse(names(V(diff_net)) %in% diff_gene_names, "red", "purple")
  
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
  plot(diff_net, layout=layout_with_fr(diff_net),
       vertex.shape=vertex_shapes,
       vertex.color=vertex_colors,
       vertex.size=vertex_sizes,
       edge.width=edge_weights,  
       vertex.label=V(diff_net)$name,
       vertex.label.cex=0.75,
       vertex.label.color="black",
       edge.arrow.size=0.5, main = "Differential network graph")
  legend("topright", legend = c("Diff Gene", "Non-Diff Gene"), pch = c(15, 15), col = c("red", "purple"))
}


# Hub genes boxplot
plot_hub_genes_boxplot <- function(expression_data_t, grouping_vector, result, dataset) {
  
  gene_strengths <- strength(result$Diff.net, weights = E(result$Diff.net)$weight)
  mean_strength <- mean(gene_strengths)
  sd_strength <- sd(gene_strengths)
  # Apply the three-sigma rule to determine hub nodes
  threshold <- mean_strength + 3 * sd_strength
  hub_genes <- names(which(gene_strengths > threshold))
  
  print(hub_genes)
  hub_gene_data <- dataset$X[, hub_genes]
  
  # Figure
  df_for_plot <- data.frame(
    Expression = as.vector(hub_gene_data),
    Gene = rep(hub_genes, each = nrow(hub_gene_data)),
    Group = rep(dataset$group, times = length(hub_genes))
  )
  
  ggplot(df_for_plot, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("green", "yellow")) +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    labs(x = "Gene", y = "Gene expression levels", title = "Hub genes boxplot") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' @Encanpsulation_internal
# handling the paper's dataset
replicate_analysis_internal <- function(dataset, condition1, condition2, output_path, pdf_name_suffix) {
  
  eval(parse(text = paste0('data("', dataset, '")')))
  dataset <- get(dataset)
  
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  #' @figure1_A: Expression data Condition 1 and Condition 2
  pdf(file.path(output_path, paste0("Expression_data_", condition1, "_", pdf_name_suffix, ".pdf")))
  plot_individual_heatmap(expression_data_t, grouping_vector, condition1)
  dev.off()
  
  pdf(file.path(output_path, paste0("Expression_data_", condition2, "_", pdf_name_suffix, ".pdf")))
  plot_individual_heatmap(expression_data_t, grouping_vector, condition2)
  dev.off()
  
  #' @figure1_B: The difference of Partial Correlation
  pdf(file.path(output_path, paste0("Difference_of_Partial_Correlation_", pdf_name_suffix, ".pdf")))
  print(plot_partial_cor_diff_matrix(results$diff.edge.weight))
  dev.off()
  
  #' @figure3: Differential network graph
  diff_gene_names <- names(results$diff.gene[results$diff.gene == 1])
  pdf(file.path(output_path, paste0("Differential_network_graph_", pdf_name_suffix, ".pdf")))
  print(plot_diff_network(results$Diff.net, diff_gene_names))
  dev.off()
  
  #' @figure4: Hub genes boxplot
  pdf(file.path(output_path, paste0("Hub_genes_boxplot_", pdf_name_suffix, ".pdf")))
  print(plot_hub_genes_boxplot(expression_data_t, grouping_vector, results, dataset))
  dev.off()
}


#' @Encanpsulation_external
# handling the reference dataset
replicate_analysis_external <- function(input1, input2, condition1, condition2, output_path, pdf_name_suffix) {
  
  ref_output <- read_and_analyze(input1, input2)
  dataset <- ref_output$data
  results <- ref_output$results
  
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  #' @figure1_A: Expression data for Condition 1 and Condition 2
  pdf(file.path(output_path, paste0("Expression_data_", condition1, "_", pdf_name_suffix, ".pdf")))
  plot_individual_heatmap(expression_data_t, grouping_vector, condition1)
  dev.off()
  
  pdf(file.path(output_path, paste0("Expression_data_", condition2, "_", pdf_name_suffix, ".pdf")))
  plot_individual_heatmap(expression_data_t, grouping_vector, condition2)
  dev.off()
  
  #' @figure1_B: The difference of Partial Correlation
  pdf(file.path(output_path, paste0("Difference_of_Partial_Correlation_", pdf_name_suffix, ".pdf")))
  plot_partial_cor_diff_matrix(results$diff.edge.weight)
  dev.off()
  
  #' @figure3: Differential network graph
  diff_gene_names <- names(results$diff.gene[results$diff.gene == 1])
  pdf(file.path(output_path, paste0("Differential_network_graph_", pdf_name_suffix, ".pdf")))
  plot_diff_network(results$Diff.net, diff_gene_names)
  dev.off()
  
  #' @figure4: Hub genes boxplot
  pdf(file.path(output_path, paste0("Hub_genes_boxplot_", pdf_name_suffix, ".pdf")))
  print(plot_hub_genes_boxplot(expression_data_t, grouping_vector, results, dataset))
  dev.off()
}


#' @test
args <- commandArgs(trailingOnly = TRUE)

operation_type <- args[1]
# chNet package dataset
if (operation_type == "internal") {
  dataset <- args[2]
  condition1 <- args[3]
  condition2 <- args[4]
  output_path <- args[5]
  pdf_name <- args[6]
  replicate_analysis_internal(dataset, condition1, condition2, output_path, pdf_name)
} else if (operation_type == "external") {
# reference dataset
  input1 <- args[2]
  input2 <- args[3]
  condition1 <- args[4]
  condition2 <- args[5]
  output_path <- args[6]
  pdf_name <- args[7]
  replicate_analysis_external(input1, input2, condition1, condition2, output_path, pdf_name)
}

## run these code terminal arguments. Change the directory accordingly
# Rscript replicate.R internal "TCGA.BRCA" "Basal" "LumA" "Output/replicate/TCGA.BRCA" TCGA.BRCA

# Rscript replicate.R internal GSE13159.AML "cancer" "normal" "output/replicate/GSE13159.AML" GSE13159.AML


# Rscript replicate.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/500_genes" gene500

# Rscript replicate.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/1000_genes" gene1000

# Rscript replicate.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/2500_genes" gene2500


## terminal path for my personal test
# cd F:/FAU/Semester\ 4/chNet_BioNet/R