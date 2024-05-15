library(chNet)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(corrplot)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(plotly)


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
  
  return(invisible(list(data = list(X = X_matrix, genes = rownames(X_matrix), group = group_vector), results = analysis_results)))
}


#' @Visualization
## Histogram of Degree Distribution 
network_histogram <- function(diff_net) {
  degree_distribution <- degree(diff_net)
  df_degrees <- data.frame(Gene = names(degree_distribution), Degree = degree_distribution)
  p <- ggplot(df_degrees, aes(x = Degree)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = "Degree Distribution of Network Nodes", x = "Number of Connections", y = "Frequency")
  return(p)  
}


## Clustering Plot
network_clustering <- function(diff_net) {
  # Calculate gene strengths and determine hub genes using the three-sigma rule
  gene_strengths <- strength(diff_net, mode = "all", weights = NULL)
  mean_strength <- mean(gene_strengths)
  sd_strength <- sd(gene_strengths)
  threshold <- mean_strength + 3 * sd_strength
  hub_genes <- names(which(gene_strengths > threshold))
  # filter the non connected
  non_isolated_nodes <- V(diff_net)[igraph::degree(diff_net, mode = "all") > 0]
  diff_net <- induced_subgraph(diff_net, non_isolated_nodes)
  # filtering nodes
  if (length(V(diff_net)) > 400) {
    degree_threshold <- quantile(igraph::degree(diff_net, mode = "all"), 0.9)
    high_degree_nodes <- V(diff_net)[igraph::degree(diff_net, mode = "all") >= degree_threshold]
    diff_net <- induced_subgraph(diff_net, high_degree_nodes)
  }
  degrees <- igraph::degree(diff_net, mode = "all")
  max_size <- 13 
  min_size <- 5  
  size_range <- max_size - min_size
  vertex_sizes <- (degrees - min(degrees)) / max(degrees) * size_range + min_size
  vertex_sizes[names(vertex_sizes) %in% hub_genes] <- max_size
  
  vertex_shapes <- ifelse(names(V(diff_net)) %in% hub_genes, "square", "circle")
  
  communities <- cluster_louvain(diff_net)
  vertex_colors <- rainbow(max(communities$membership))[communities$membership]
  
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
       edge.arrow.size=0.5, main = "Clustering of Differential network graph")
}


## Boxplot for the most significant genes
hub_gene_boxplot <- function(dataset, hub_gene_data, top_hub_genes){
  df_for_plot <- data.frame(
    Expression = as.vector(hub_gene_data),
    Gene = rep(top_hub_genes, each = nrow(hub_gene_data)),
    Group = rep(dataset$group, times = length(top_hub_genes))
  )
  
  ggplot(df_for_plot, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("blue", "red")) +
    theme_minimal() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Expression Levels", title = "Hub Gene Expression Boxplot")
}


## Violin for the most significant genes
hub_gene_violin <- function(dataset, hub_gene_data, top_hub_genes){
  df_for_plot <- data.frame(
    Expression = as.vector(hub_gene_data),
    Gene = rep(top_hub_genes, each = nrow(hub_gene_data)),
    Group = rep(dataset$group, times = length(top_hub_genes))
  )
  
  ggplot(df_for_plot, aes(x = Gene, y = Expression, fill = Group)) +
    geom_violin() +
    theme_minimal() +
    labs(x = "Gene", y = "Expression Level", title = "Hub Gene Expression Violin Plot") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


## Heatmap for the most significant genes
hub_gene_heatmap <- function(dataset, top_hub_genes, group_vector) {
  
  valid_genes_mask <- colnames(dataset$X) %in% top_hub_genes
  valid_top_hub_genes <- colnames(dataset$X)[valid_genes_mask]
  
  hub_gene_data <- dataset$X[, valid_top_hub_genes, drop = FALSE]
  z_scores <- t(scale(t(as.matrix(hub_gene_data))))
  
  # Handle NA/Inf values
  z_scores[is.na(z_scores)] <- 0
  z_scores[is.infinite(z_scores)] <- 0
  
  gene_variances <- apply(z_scores, 1, var)
  if (length(gene_variances) > 600) {
    most_variable_genes <- names(sort(gene_variances, decreasing = TRUE)[1:600])
    z_scores <- z_scores[most_variable_genes,]
  }
  
  relevant_sample_indices <- which(colnames(dataset$X) %in% colnames(z_scores))
  adjusted_group_vector <- group_vector[relevant_sample_indices]
  
  group_colors <- rainbow(length(unique(na.omit(adjusted_group_vector))))
  names(group_colors) <- unique(na.omit(adjusted_group_vector))
  if (any(is.na(adjusted_group_vector))) {
    group_colors["NA"] <- "grey"  # Explicitly assign grey to NA
  }

  group_annotation <- data.frame(Group = adjusted_group_vector)
  ha <- HeatmapAnnotation(df = group_annotation, col = list(Group = group_colors))
  
  Heatmap(z_scores, name = "Z-score",
          show_row_names = FALSE, show_column_names = TRUE,
          cluster_rows = TRUE, cluster_columns = TRUE,
          row_title = "Top Hub Genes Expression", column_title = "Samples",
          top_annotation = ha,
          color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
}


## Correlation for the most significant genes
hub_gene_correlation <- function(dataset, hub_gene_data, top_hub_genes){
  cor_matrix <- cor(hub_gene_data)
  corrplot(cor_matrix, method = "circle", type = "upper",
           order = "hclust",
           tl.col = "black",    
           tl.srt = 45,       
           main = "Hub Gene Expression Correlation Plot", mar = c(4, 4, 2, 2))
}


# Gene Ontology Enrichment - Molecular Functions of hub genes  
Gene_Ontology_Enrichment_MF <- function(dataset, numHubGenes) {
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  diff_net <- results$Diff.net
  gene_strengths <- strength(diff_net, weights = E(diff_net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:numHubGenes])
  
  entrez_ids <- mapIds(org.Hs.eg.db, keys = top_hub_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  if (length(entrez_ids) == 0) {
    print("No Entrez IDs found for top hub genes.")
    return(NULL)
  }
  
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
      
      p <- ggplot(top_df, aes(x = reorder(label, pvalue), y = Count)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = Count), vjust = -0.3, color = "black", size = 3.5) +
        labs(x = "Molecular Functions and Genes", y = "Number of Connected Genes", title = "Top Hub Genes GO Enrichment Analysis - Molecular Functions") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 12))
      
      return(p)
    } else {
      print("No significant GO terms found for MF ontology.")
      return(NULL)
    }
  } else {
    print("No significant GO terms found for MF ontology.")
    return(NULL)
  }
}

# Gene Ontology Enrichment - Biological Process of hub genes  
Gene_Ontology_Enrichment_BP <- function(dataset, numHubGenes) {
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R = 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  diff_net <- results$Diff.net
  gene_strengths <- strength(diff_net, weights = E(diff_net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:numHubGenes])
  
  entrez_ids <- mapIds(org.Hs.eg.db, keys = top_hub_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  if (length(entrez_ids) == 0) {
    print("No Entrez IDs found for top hub genes.")
    return(NULL)
  }
  
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
      
      p <- ggplot(top_df, aes(x = reorder(label, pvalue), y = Count)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = Count), vjust = -0.3, color = "black", size = 3.5) +
        labs(x = "Biological Process and Genes", y = "Number of Connected Genes", title = "Top Hub Genes GO Enrichment Analysis - Biological Process") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 12))
      
      return(p)
    } else {
      print("No significant GO terms found for BP ontology.")
      return(NULL)
    }
  } else {
    print("No significant GO terms found for BP ontology.")
    return(NULL)
  }
}

#' @Encanpsulation_internal
# handling the paper's dataset
downstream_analysis_internal <- function(dataset, condition1, condition2, output_path, pdf_name_suffix) {
  
  eval(parse(text = paste0('data("', dataset, '")')))
  dataset <- get(dataset)
  
  results <- chNet(dataset$X, dataset$group, subsampling = FALSE, R= 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  diff_net <- results$Diff.net
  
  gene_strengths <- strength(results$Diff.net, weights = E(results$Diff.net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:10])
  hub_gene_data <- dataset$X[, top_hub_genes]
  
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_BP_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_BP(dataset, numHubGenes = 10))
  dev.off()
  
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_MF_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_MF(dataset, numHubGenes = 10))
  dev.off()
  
  pdf(file.path(output_path, paste0("network_histogram_", pdf_name_suffix, ".pdf")))
  print(network_histogram(diff_net))
  dev.off()
  
  pdf(file.path(output_path, paste0("network_clustering_", pdf_name_suffix, ".pdf")))
  result <- network_clustering(diff_net)
  if (!is.null(result$plot)) { # avoid unnecessary print
    print(result$plot)
  }
  dev.off()
  
  pdf(file.path(output_path, paste0("hub_gene_boxplot_", pdf_name_suffix, ".pdf")))
  print(hub_gene_boxplot(dataset, hub_gene_data, top_hub_genes))
  dev.off()
  
  pdf(file.path(output_path, paste0("hub_gene_correlation_", pdf_name_suffix, ".pdf")))
  hub_gene_correlation(dataset, hub_gene_data, top_hub_genes)
  dev.off()
  
  
  pdf(file.path(output_path, paste0("hub_gene_heatmap_", pdf_name_suffix, ".pdf")))
  print(hub_gene_heatmap(dataset, top_hub_genes, grouping_vector))
}


#' @Encanpsulation_external
# handling the reference dataset
downstream_analysis_external <- function(input1, input2, condition1, condition2, output_path, pdf_name_suffix) {
  
  ref_output <-  read_and_analyze(input1, input2)
  dataset <- ref_output$data
  results <- ref_output$results
  
  expression_data_t <- t(dataset$X)
  grouping_vector <- dataset$group
  
  diff_net <- results$Diff.net
  
  gene_strengths <- strength(results$Diff.net, weights = E(results$Diff.net)$weight)
  top_hub_genes <- names(sort(gene_strengths, decreasing = TRUE)[1:10])
  hub_gene_data <- dataset$X[, top_hub_genes]
  
  pdf(file.path(output_path, paste0("network_histogram_", pdf_name_suffix, ".pdf")))
  print(network_histogram(diff_net))
  dev.off()
  
  pdf(file.path(output_path, paste0("network_clustering_", pdf_name_suffix, ".pdf")))
  result <- network_clustering(diff_net)
  if (!is.null(result$plot)) { # avoid unnecessary print
    print(result$plot)
  }
  dev.off()
  
  pdf(file.path(output_path, paste0("hub_gene_boxplot_", pdf_name_suffix, ".pdf")))
  print(hub_gene_boxplot(dataset, hub_gene_data, top_hub_genes))
  dev.off()
  
  pdf(file.path(output_path, paste0("hub_gene_correlation_", pdf_name_suffix, ".pdf")))
  hub_gene_correlation(dataset, hub_gene_data, top_hub_genes)
  dev.off()
  
  pdf(file.path(output_path, paste0("hub_gene_heatmap_", pdf_name_suffix, ".pdf")))
  print(hub_gene_heatmap(dataset, top_hub_genes, grouping_vector))
  dev.off()
  
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_BP_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_BP(dataset, numHubGenes = 10))
  dev.off()
  
  pdf(file.path(output_path, paste0("gene_ontology_enrichment_MF_", pdf_name_suffix, ".pdf")), width = 11, height = 10)
  print(Gene_Ontology_Enrichment_MF(dataset, numHubGenes = 10))
  dev.off()
}


#' @test
args <- commandArgs(trailingOnly = TRUE)

operation_type <- args[1]
if (operation_type == "internal") {
  dataset <- args[2]
  condition1 <- args[3]
  condition2 <- args[4]
  output_path <- args[5]
  pdf_name <- args[6]
  downstream_analysis_internal(dataset, condition1, condition2, output_path, pdf_name)
} else if (operation_type == "external") {
  input1 <- args[2]
  input2 <- args[3]
  condition1 <- args[4]
  condition2 <- args[5]
  output_path <- args[6]
  pdf_name <- args[7]
  downstream_analysis_external(input1, input2, condition1, condition2, output_path, pdf_name)
}

## run these code terminal arguments. Change the directory accordingly
# Rscript downstream.R internal "TCGA.BRCA" "Basal" "LumA" "output/downstream/TCGA.BRCA" TCGA.BRCA

# Rscript downstream.R internal GSE13159.AML "cancer" "normal" "output/downstream/GSE13159.AML" GSE13159.AML


# Rscript downstream.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/500_genes" gene500


# Rscript downstream.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/1000_genes" gene1000

# Rscript downstream.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/2500_genes" gene2500


## terminal path for my personal test
# cd F:/FAU/Semester\ 4/chNet_BioNet/R