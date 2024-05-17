# Installing the DiffCorr Package
install.packages("BiocManager")
BiocManager::install("http://bioconductor.org/biocLite.R", force = TRUE)
BiocManager::install(c("pcaMethods", "multtest"), force = TRUE)
BiocManager::install("GEOquery", force = TRUE)
BiocManager::install("affy", force = TRUE)
BiocManager::install("genefilter", force = TRUE)
BiocManager::install("GOstats", force = TRUE)
BiocManager::install("ath1121501.db", force = TRUE)
install.packages("spatstat", force = TRUE)
install.packages("igraph", force = TRUE)
library(DiffCorr)


## Reading the Golub dataset
data(golub)
dim(golub)

## Clusters on each subset
hc.mol1 <- cluster.molecule(golub[, 1:27], "pearson", "average")  ## ALL (27 samples)
hc.mol2 <- cluster.molecule(golub[, 28:38], "pearson", "average") ## AML (11 samples)


## Cut the tree at a correlation of 0.6 using cutree function
g1 <- cutree(hc.mol1, h=0.4)
g2 <- cutree(hc.mol2, h=0.4)


##
res1 <- get.eigen.molecule(data = golub, groups = g1)
res2 <- get.eigen.molecule(data = golub, groups = g2)

###################################
## Visualizing module networks
###################################
gg1 <- get.eigen.molecule.graph(res1)
plot(gg1, layout=layout.fruchterman.reingold(gg1), main = "ALL")
write.modules(g1, res1, outfile="module1_list.txt")

gg2 <- get.eigen.molecule.graph(res2)
plot(gg2, layout=layout.fruchterman.reingold(gg2), main = "AML")
write.modules(g2, res2, outfile="module2_list.txt")


###################################
## You can examine the relationship between modules
###################################
for (i in 1:length(res1$eigen.molecules)) {
  for (j in 1: length(res2$eigen.molecules)) {
    r <- cor(res1$eigen.molecules[[i]],res2$eigen.molecules[[j]], method="spearman")
    if (abs(r) > 0.8) {
      print(paste("(i, j): ", i, " ", j, sep=""))
      print(r)
    }
  }
}

cor(res1$eigen.molecules[[2]],res2$eigen.molecules[[8]], method="spearman")
plot(res1$eigen.molecules[[2]], res2$eigen.molecules[[8]], main = "ALL")
plot(res1$eigen.molecules[[21]], res2$eigen.molecules[[24]], main = "AML")

####################################
## Examine groups of interest graphically
## look at groups 21 and 24 
####################################
plotDiffCorrGroup(golub, g1, g2, 21, 24, 1:27, 28:38,
                  scale.center=TRUE, scale.scale=TRUE,
                  ylim=c(-5,5))

####################################
## Export the results (FDR < 0.05)
####################################
comp.2.cc.fdr(output.file="Result.txt", golub[,1:27], golub[,28:38], threshold=0.05)

