library(chNet)
library(igraph)


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


save_differential_network <- function(result, output_dir, file_name) {
  if (!inherits(result$Diff.net, "igraph")) {
    stop("Diff.net is not an igraph object.")
  }
  
  library(igraph) 
  
  diff_net_edges <- as.data.frame(as_edgelist(result$Diff.net, names = TRUE))
  
  diff_net_edges$weight <- E(result$Diff.net)$weight
  
  diff_status <- result$diff.gene
  
  label_edge_condition <- function(edge_from, edge_to, diff_status) {
    condition <- ifelse(diff_status[edge_from] == 1 && diff_status[edge_to] == 1, "differential", "non-differential")
    return(condition)
  }
  
  diff_net_edges$condition <- mapply(label_edge_condition, 
                                     diff_net_edges$V1, 
                                     diff_net_edges$V2, 
                                     MoreArgs = list(diff_status = diff_status))
  
  diff_net_edges <- diff_net_edges[, c("V1", "V2", "condition", "weight")]
  names(diff_net_edges) <- c("regulator", "target", "condition", "weight")
  
  
  output_file_path <- file.path(output_dir, file_name)
  write.table(diff_net_edges, output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(paste("Network file saved to", output_file_path))
}


#' @Encanpsulation_internal
# handling the paper's dataset
run_chNet_analysis_internal <- function(data_set_name, output_dir, file_name) {
  
  eval(parse(text = paste0('data("', data_set_name, '")')))
  dataset <- get(data_set_name)
  
  result <- chNet(
    X = dataset$X, 
    group = dataset$group, 
    subsampling = FALSE, R= 20, lambar = 2.825, parallel = FALSE, nCpus = 4)
  
  save_differential_network(result, output_dir, file_name)
}


#' @Encanpsulation_external
# handling the reference dataset
run_chNet_analysis_external <- function(input1, input2, output_dir, file_name) {
  
  ref_output <- read_and_analyze(input1, input2)
  ref_dataset <- ref_output$combined_data
  ref_result <- ref_output$results
  
  save_differential_network(ref_result, output_dir, file_name)
}


#' @test
args <- commandArgs(trailingOnly = TRUE)

operation_type <- args[1]
# chNet package dataset
if (operation_type == "internal") {
  data_set_name <- args[2]
  output_dir <- args[3]
  file_name <- args[4]
  run_chNet_analysis_internal(data_set_name, output_dir, file_name)
} else if (operation_type == "external") {
# reference dataset
  input1 <- args[2]
  input2 <- args[3]
  output_dir <- args[4]
  file_name <- args[5]
  run_chNet_analysis_external(input1, input2, output_dir, file_name)
}

