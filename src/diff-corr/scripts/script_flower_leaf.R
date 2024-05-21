# Load necessary packages
library(affy)
library(genefilter)
library(igraph)
library(spatstat)
library(cluster)
library(DiffCorr)

# Create a logging function
log_message <- function(message, log_file = "diffcorr_analysis_leaf.log") {
  timestamp <- Sys.time()
  formatted_message <- paste0("[", timestamp, "] ", message, "\n")
  cat(formatted_message, file = log_file, append = TRUE)
}

# Log the start of the process
log_message("Starting DiffCorr analysis for LEaf dataset")

# Downloading the Transcriptome Data set
log_message("Downloading GSE5632 data")
data <- getGEOSuppFiles("GSE5632") 
untar("GSE5632/GSE5632_RAW.tar", exdir="GSE5632")

log_message("Downloading GSE5630 data")
data <- getGEOSuppFiles("GSE5630") 
untar("GSE5630/GSE5630_RAW.tar", exdir="GSE5630")

# Filtering - RMA
log_message("Performing RMA normalization on GSE5630 data")
tgt <- list.files("./GSE5630", pattern="*.CEL.gz", full.names=TRUE)
eset.GSE5630 <- justRMA(filenames=tgt)

log_message("Performing RMA normalization on GSE5632 data")
tgt2 <- list.files("./GSE5632", pattern="*.CEL.gz", full.names=TRUE)
eset.GSE5632 <- justRMA(filenames=tgt2) 

# Discard all control probes
log_message("Discarding control probes")
rmv <- c(grep("AFFX", rownames(eset.GSE5632)), grep("s_at", rownames(eset.GSE5632)), grep("x_at", rownames(eset.GSE5632)))
eset.GSE5632 <- eset.GSE5632[-rmv,] 
eset.GSE5630 <- eset.GSE5630[-rmv,] 

# Filter function for the expression level and the coefficient of variation
log_message("Applying filter function for expression level and coefficient of variation")
e.mat <- 2 ^ exprs(eset.GSE5632)
ffun <- filterfun(pOverA(0.2, 100), cv(0.5, 10))
filtered <- genefilter(e.mat, ffun)
eset.GSE5632.sub <- log2(e.mat[filtered, ])

e.mat <- 2 ^ exprs(eset.GSE5630)
filtered <- genefilter(e.mat, ffun)
eset.GSE5630.sub <- log2(e.mat[filtered, ])

# Identify common probe sets between the two data sets
log_message("Identifying common probe sets between GSE5630 and GSE5632")
comm <- intersect(rownames(eset.GSE5632.sub), rownames(eset.GSE5630.sub))
eset.GSE5632.sub <- eset.GSE5632.sub[comm, ]
eset.GSE5630.sub <- eset.GSE5630.sub[comm, ]

data1 <- cbind(eset.GSE5632.sub, eset.GSE5630.sub)

# Clustering on each subset
log_message("Clustering on each subset")
hc.flowers <- cluster.molecule(data1[, 1:66], method = "pearson", linkage = "average")
hc.leaves <- cluster.molecule(data1[, 67:126], method = "pearson", linkage = "average")

# Cut the tree at a correlation of 0.6 using cutree function
log_message("Cutting tree at correlation of 0.4")
g1 <- cutree(hc.flowers, h = 0.4)
g2 <- cutree(hc.leaves, h = 0.4)

log_message("Getting eigen molecules")
res1 <- get.eigen.molecule(data1, groups = g1, whichgroups = c(1:10), methods = "svd", n = 2)
res2 <- get.eigen.molecule(data1, groups = g2, whichgroups = c(11:20), methods = "svd", n = 2)

# Visualizing module networks
log_message("Visualizing module networks for GSE530")
gg1 <- get.eigen.molecule.graph(res1)
plot(gg1, layout = layout.fruchterman.reingold(gg1), main = "GSE530")
write.modules(g1, res1, outfile = "module1_list.txt")

log_message("Visualizing module networks for GSE532")
gg2 <- get.eigen.molecule.graph(res2)
plot(gg2, layout = layout.fruchterman.reingold(gg2), main = "GSE532")
write.modules(g2, res2, outfile = "module2_list.txt")

log_message("Plotting differential correlation group")
plotDiffCorrGroup(data1, g1, g2, 1, 11, 1:66, 67:126, scale.center = TRUE, scale.scale = TRUE, ylim = c(-5, 5))

# Export the results
log_message("Exporting results with FDR < 0.05")
comp.2.cc.fdr(output.file = "Transcript_DiffCorr_res.txt", data1[, 1:66], data1[, 67:126], threshold = 0.05)

log_message("DiffCorr analysis completed")
