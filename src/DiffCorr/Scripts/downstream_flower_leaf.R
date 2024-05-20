# Installing the DiffCorr Package
install.packages("BiocManager")
BiocManager::install("http://bioconductor.org/biocLite.R")
BiocManager::install(c("pcaMethods", "multtest"), force = TRUE)
library(DiffCorr)
BiocManager::install("GEOquery", force = TRUE)
BiocManager::install("affy", force = TRUE)
BiocManager::install("genefilter", force = TRUE)
BiocManager::install("GOstats", force = TRUE)
BiocManager::install("ath1121501.db", force = TRUE)
install.packages("spatstat", force = TRUE)
install.packages("igraph", force = TRUE)

# Downloading the Transcriptome Data set
library("GEOquery")
data <- getGEOSuppFiles("GSE5632") 
untar("GSE5632/GSE5632_RAW.tar", exdir="GSE5632")
data <- getGEOSuppFiles("GSE5630") 
untar("GSE5630/GSE5630_RAW.tar", exdir="GSE5630")

# Filtering - RMA
library(affy)
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
library(genefilter)
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

# correlation matrices for each data set by Spearman’s rank- order correlation
GSE5632.cor <- cor(t(eset.GSE5632.sub), method="spearman") 
GSE5630.cor <- cor(t(eset.GSE5630.sub), method="spearman")

# Visualization on a pseudo-color heatmap
library(spatstat)
par(mfrow=c(1,2)) 
plot(im(GSE5632.cor[nrow(GSE5632.cor):1,]),col=cm.colors(256), main="GSE5632") 
plot(im(GSE5630.cor[nrow(GSE5630.cor):1,]),col=cm.colors(256), main="GSE5630")

# Construction of the co-expression network
library(igraph)
g1 <- graph.adjacency(GSE5632.cor, weighted=TRUE, mode="lower")
g1 <- delete.edges(g1, E(g1)[ weight < 0.95 ])
g1 <- igraph::simplify(g1, remove.multiple = TRUE, remove.loops = TRUE)
g1 <- delete.vertices(g1, which(igraph::degree(g1)<1)) 
plot(g1, vertex.size=3, edge.width=3, vertex.color=ifelse(igraph::degree(g1)>20,"Magenta","Green"), vertex.label="", layout=layout.kamada.kawai, main="GSE5632")

g2 <- graph.adjacency(GSE5630.cor, weighted=TRUE, mode="lower")
g2 <- delete.edges(g2, E(g2)[ weight < 0.95 ])
g2 <- igraph::simplify(g2, remove.multiple = TRUE,remove.loops = TRUE)
g2 <- delete.vertices(g2, which(igraph::degree(g2)<1)) 
plot(g2, vertex.size=3, edge.width=3, vertex.color=ifelse(igraph::degree(g2)>20,"Magenta","Green"), vertex.label="", layout=layout.kamada.kawai, main="GSE5630")

# Cytoscape
write.graph(g1, "g1forcy.gml", format="gml")
write.graph(g2, "g2forcy.gml", format="gml")

# Graph Clustering
g1.fc <- fastgreedy.community(g1)
sizes(g1.fc)

g2.fc <- fastgreedy.community(g2)
sizes(g2.fc)

# access each module member
mod1 <- membership(g1.fc)[membership(g1.fc)==1]
mod1.p <- names(mod1)

mod2 <- membership(g1.fc)[membership(g1.fc)==2] 
mod2.p <- names(mod2)

mod3 <- membership(g1.fc)[membership(g1.fc)==3] 
mod3.p <- names(mod3)

# Gene Ontology (GO) term enrichment analyses were performed to assess cluster fidelity
# Gene Ontology Enrichment Analysis
library(GOstats)
library(GO.db)
library(ath1121501.db)
ls("package:ath1121501.db")
x <- ath1121501ACCNUM
mapped.probes <- mappedkeys(x)
length(mapped.probes)
geneUniv <- AnnotationDbi::as.list(x[mapped.probes])
mod1.p.gene <- unique(unlist(AnnotationDbi::as.list(x[mod1.p])))
mod2.p.gene <- unique(unlist(AnnotationDbi::as.list(x[mod2.p])))
mod3.p.gene <- unique(unlist(AnnotationDbi::as.list(x[mod3.p])))
hgCutoff <- 0.0001

params <- new("GOHyperGParams",geneIds=mod1.p.gene, universeGeneIds=geneUniv, annotation="ath1121501", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
hgOver <- hyperGTest(params)
df <- summary(hgOver) 
names(df)
pvalues(hgOver)[1:3]
htmlReport(hgOver, file="res_mod1.html")

params <- new("GOHyperGParams", geneIds=mod2.p.gene, universeGeneIds=geneUniv, annotation="ath1121501", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
hgOver <- hyperGTest(params)
htmlReport(hgOver, file="/Users/rameshsubramani/Downloads/DiffCorr/GSE/res_mod2.html")

params <- new("GOHyperGParams", geneIds=mod3.p.gene, universeGeneIds=geneUniv, annotation="ath1121501", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
hgOver <- hyperGTest(params)
htmlReport(hgOver, file="/Users/rameshsubramani/Downloads/DiffCorr/GSE/res_mod3.html")


# Differential Correlation Analysis by DiffCorr Package
# Calculation of Differential Co-Expression between Organs in Arabidopsis
dim(eset.GSE5632.sub)
dim(eset.GSE5630.sub)
data1 <- cbind(eset.GSE5632.sub, eset.GSE5630.sub) 
?cluster.molecule
hc.flowers <- cluster.molecule(data1[, 1:66],method="pearson", linkage="average")
hc.leaves <- cluster.molecule(data1[, 67:126],method="pearson", linkage="average")  
## Cut the tree at a correlation of 0.6 using cutree function
#library(dynamicTreeCut)
g1 <- cutree(hc.flowers, h=0.4) 
g2 <- cutree(hc.leaves, h=0.4)
res1 <- get.eigen.molecule(data1, groups=g1,whichgroups=c(1:10), methods="svd", n=2)
res2 <- get.eigen.molecule(data1, groups=g2,whichgroups=c(11:20), methods="svd", n=2)
gg1 <- get.eigen.molecule.graph(res1)
plot(gg1, layout=layout.fruchterman.reingold(gg1), main="GSE5630")
write.modules(g1, res1, outfile="module1_list.txt") 

gg2 <- get.eigen.molecule.graph(res2)
plot(gg2, layout=layout.fruchterman.reingold(gg2), main="GSE5632") 
write.modules(g2, res2, outfile="module2_list.txt")

write.graph(gg1, "tmp1.ncol", format="ncol") 
write.graph(gg2, "tmp2.ncol", format="ncol") 
tmp1 <- read.table("tmp1.ncol")
tmp2 <- read.table("tmp2.ncol")
tmp1$V3 <- "pp" 
tmp2$V3 <- "pp"
tmp1$V4 <- "gg1forcy"
tmp2$V4 <- "gg2forcy"
tmp1 <- tmp1[, c("V4", "V1", "V3", "V2")]
tmp2 <- tmp2[, c("V4", "V1", "V3", "V2")]
write.table(tmp1, file="gg1forcy.nnf", row.names=FALSE,col.names=FALSE)
write.table(tmp2, file="gg2forcy.nnf", row.names=FALSE,col.names=FALSE)
module1_list <- read.table("module1_list.txt", skip=1) 
module2_list <- read.table("module2_list.txt", skip=1) 
module1_list$V1 <- sub("̂", "Module", module1_list$V1)
module2_list$V1 <- sub("̂", "Module", module2_list$V1 -10)
write.table(module1_list, file="gg1forcy.nnf", append=TRUE,row.names=FALSE, col.names=FALSE)
write.table(module2_list, file="gg2forcy.nnf", append=TRUE, row.names=FALSE, col.names=FALSE)