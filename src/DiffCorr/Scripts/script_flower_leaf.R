# Installing the DiffCorr Package
install.packages("BiocManager")
BiocManager::install("http://bioconductor.org/biocLite.R")
BiocManager::install(c("pcaMethods", "multtest"), force = TRUE)
BiocManager::install("GEOquery", force = TRUE)
BiocManager::install("affy", force = TRUE)
BiocManager::install("genefilter", force = TRUE)
BiocManager::install("GOstats", force = TRUE)
BiocManager::install("ath1121501.db", force = TRUE)
install.packages("spatstat", force = TRUE)
install.packages("igraph", force = TRUE)

library(DiffCorr)
library("GEOquery")
library(affy)
library(genefilter)

# Downloading the Transcriptome Data set
data <- getGEOSuppFiles("GSE5632") 
untar("GSE5632/GSE5632_RAW.tar", exdir="GSE5632")
data <- getGEOSuppFiles("GSE5630") 
untar("GSE5630/GSE5630_RAW.tar", exdir="GSE5630")

# Filtering - RMA
tgt <- list.files("./GSE5630", pattern="*.CEL.gz", full.names=TRUE)
eset.GSE5630 <- justRMA(filenames=tgt)

tgt2 <- list.files("./GSE5632", pattern="*.CEL.gz", full.names=TRUE)
eset.GSE5632 <- justRMA(filenames=tgt2) 

# discard all control probes
rmv <- c(grep("AFFX", rownames(eset.GSE5632)), grep("s_at", rownames(eset.GSE5632)),grep("x_at", rownames(eset.GSE5632)))
eset.GSE5632 <- eset.GSE5632[-rmv,] 
eset.GSE5630 <- eset.GSE5630[-rmv,] 

# We use a filter function for the expression level and the coefficient of variation. 
# The ratio of the standard deviation and the mean of a gene’s expression values across all samples must be higher than a given threshold.
e.mat <- 2 ^ exprs(eset.GSE5632)
ffun <- filterfun(pOverA(0.2, 100), cv(0.5, 10))
filtered <- genefilter(e.mat,ffun)
eset.GSE5632.sub <- log2(e.mat[filtered, ])

e.mat <- 2 ^ exprs(eset.GSE5630)
ffun <- filterfun(pOverA(0.2, 100), cv(0.5, 10))
filtered <- genefilter(e.mat,ffun)
eset.GSE5630.sub <- log2(e.mat[filtered, ])

# Identify common probe sets between the two data sets.
comm <- intersect(rownames(eset.GSE5632.sub), rownames(eset.GSE5630.sub))
eset.GSE5632.sub <- eset.GSE5632.sub[comm, ]
eset.GSE5630.sub <- eset.GSE5630.sub[comm, ]

data1 <- cbind(eset.GSE5632.sub, eset.GSE5630.sub) 
hc.flowers <- cluster.molecule(data1[, 1:66],method="pearson", linkage="average")
hc.leaves <- cluster.molecule(data1[, 67:126],method="pearson", linkage="average")  
## Cut the tree at a correlation of 0.6 using cutree function
#library(dynamicTreeCut)
g1 <- cutree(hc.flowers, h=0.4) 
g2 <- cutree(hc.leaves, h=0.4)
res1 <- get.eigen.molecule(data1, groups=g1, whichgroups=c(1:10), methods="svd", n=2)
res2 <- get.eigen.molecule(data1, groups=g2,whichgroups=c(11:20), methods="svd", n=2)
gg1 <- get.eigen.molecule.graph(res1)
plot(gg1, layout=layout.fruchterman.reingold(gg1), main="GSE530")
write.modules(g1, res1, outfile="module1_list.txt") 

gg2 <- get.eigen.molecule.graph(res2)
plot(gg2, layout=layout.fruchterman.reingold(gg2), main="GSE532") 
write.modules(g2, res2, outfile="module2_list.txt")

plotDiffCorrGroup(data1, g1, g2, 1, 11, 1:66, 67:126,scale.center=TRUE, scale.scale=TRUE,ylim=c(-5,5))

comp.2.cc.fdr(output.file="Transcript_DiffCorr_res.txt",data1[,1:66], data1[,67:126], threshold=0.05)
