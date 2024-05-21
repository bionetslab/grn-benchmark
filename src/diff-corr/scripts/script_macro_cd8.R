# Load necessary packages
library(affy)
library(genefilter)
library(igraph)
library(spatstat)
library(cluster)
library(DiffCorr)

# Create a logging function
log_message <- function(message, log_file = "diffcorr_analysis_macro.log") {
  timestamp <- Sys.time()
  formatted_message <- paste0("[", timestamp, "] ", message, "\n")
  cat(formatted_message, file = log_file, append = TRUE)
}

# Log the start of the process
log_message("Starting DiffCorr analysis for Macro")

# Set the working directory to the script's directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define the function
process_datasets <- function(input_file1, input_file2, output_file) {
  log_message("Reading the data")
  
  # Read the data into R
  eset_macro <- read.table(input_file1, header = TRUE, sep = "\t")
  row.names(eset_macro) <- eset_macro[, 1]
  eset_macro <- eset_macro[, -1]
  
  eset_cd8 <- read.table(input_file2, header = TRUE, sep = "\t")
  row.names(eset_cd8) <- eset_cd8[, 1]
  eset_cd8 <- eset_cd8[, -1]
  
  data1 <- cbind(eset_macro, eset_macro) 
  dim(data1)
  
  log_message("Clustering on each subset")
  
  hc.macro <- cluster.molecule(data1[1:100], method = "pearson", linkage = "average")  
  hc.cd8 <- cluster.molecule(data1[101:200], method = "pearson", linkage = "average")
  gg.macro <- cutree(hc.macro, h = 0.4) 
  gg.cd8 <- cutree(hc.cd8, h = 0.4)
  table(gg.macro[table(gg.macro) != 1])
  table(gg.macro)
  table(gg.cd8)
  
  log_message("Getting eigen molecules")
  
  res.macro <- get.eigen.molecule(data1, groups = gg.macro, methods = "svd", n = 2)
  res.cd8 <- get.eigen.molecule(data1, groups = gg.cd8, methods = "svd", n = 2)
  
  log_message("Visualizing module networks")
  
  macro_g <- get.eigen.molecule.graph(res.macro)
  plot(macro_g, layout = layout.fruchterman.reingold(macro_g), main = "Macro")
  # write.modules(gg.macro, res.macro, outfile = "module1_list.txt")
  
  cd8_g <- get.eigen.molecule.graph(res.cd8)
  plot(cd8_g, layout = layout.fruchterman.reingold(cd8_g), main = "CD8") 
  # write.modules(gg.cd8, res.cd8, outfile = "module2_list.txt")
  
  data2 <- cbind(eset_macro, eset_cd8)
  
  log_message("Exporting the results (FDR < 0.05)")
  
  # plotDiffCorrGroup(data2, gg.macro, gg.cd8, 1, 4, 1:5, 6:7, scale.center = TRUE, scale.scale = TRUE, ylim = c(-5, 5))
  comp.2.cc.fdr(output.file = output_file, data2[, 1:100], data2[, 101:200], threshold = 0.05)
  
  log_message("DiffCorr analysis completed")
}

# Log usage of the function
log_message("Calling process_datasets function")

tryCatch(
  {
    process_datasets(
      "out_Macrophages.tsv", 
      "out_CD8_exhausted.tsv", 
      "DiffCorr_res.txt"
    )
    log_message("process_datasets function completed successfully")
  },
  error = function(e) {
    log_message(paste("Error in process_datasets function: ", e$message))
  }
)
