# Load necessary packages
library(affy)
library(genefilter)
library(igraph)
library(spatstat)
library(cluster)
library(DiffCorr)

# Create a logging function
log_message <- function(message, log_file = "diffcorr_analysis.log") {
  timestamp <- Sys.time()
  formatted_message <- paste0("[", timestamp, "] ", message, "\n")
  cat(formatted_message, file = log_file, append = TRUE)
}

# Log the start of the process
log_message("Starting DiffCorr analysis")

# Define the function
process_datasets <- function(output_file) {
  log_message("Reading the Golub dataset")
  data(golub)
  dim(golub)
  
  log_message("Clustering on each subset")
  hc.mol1 <- cluster.molecule(golub[, 1:27], "pearson", "average")  # ALL (27 samples)
  hc.mol2 <- cluster.molecule(golub[, 28:38], "pearson", "average") # AML (11 samples)
  
  log_message("Cutting the tree at a correlation of 0.6")
  g1 <- cutree(hc.mol1, h=0.4)
  g2 <- cutree(hc.mol2, h=0.4)
  
  log_message("Getting eigen molecules")
  res1 <- get.eigen.molecule(data = golub, groups = g1)
  res2 <- get.eigen.molecule(data = golub, groups = g2)
  
  log_message("Visualizing module networks")
  gg1 <- get.eigen.molecule.graph(res1)
  plot(gg1, layout=layout.fruchterman.reingold(gg1), main = "ALL")
  # write.modules(g1, res1, outfile="module1_list.txt")
  
  gg2 <- get.eigen.molecule.graph(res2)
  plot(gg2, layout=layout.fruchterman.reingold(gg2), main = "AML")
  # write.modules(g2, res2, outfile="module2_list.txt")
  
  log_message("Examining the relationship between modules")
  for (i in 1:length(res1$eigen.molecules)) {
    for (j in 1:length(res2$eigen.molecules)) {
      r <- cor(res1$eigen.molecules[[i]], res2$eigen.molecules[[j]], method="spearman")
      if (abs(r) > 0.8) {
        log_message(paste("(i, j): ", i, " ", j, " Correlation: ", r, sep=""))
      }
    }
  }
  
  cor(res1$eigen.molecules[[2]], res2$eigen.molecules[[8]], method="spearman")
  plot(res1$eigen.molecules[[2]], res2$eigen.molecules[[8]], main = "ALL")
  plot(res1$eigen.molecules[[21]], res2$eigen.molecules[[24]], main = "AML")
  
  log_message("Exporting the results (FDR < 0.05)")
  comp.2.cc.fdr(output.file=output_file, golub[,1:27], golub[,28:38], threshold=0.05)
  
  log_message("DiffCorr analysis completed")
}

# Log usage of the function
log_message("Calling process_datasets function")

tryCatch(
  {
    process_datasets("ALL_AML_DiffCorr_res.txt")
    log_message("process_datasets function completed successfully")
  },
  error = function(e) {
    log_message(paste("Error in process_datasets function: ", e$message))
  }
)
