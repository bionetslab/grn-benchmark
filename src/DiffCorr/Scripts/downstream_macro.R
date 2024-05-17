# Load required libraries
install.packages("BiocManager")
library(affy)
library(genefilter)
library(igraph)
library(spatstat)
library(cluster)
library(DiffCorr)

# Read the data into R
eset_macro <- read.table("/Users/rameshsubramani/Downloads/DiffCorr/Macro_CD8/500_genes/out_Macrophages.tsv", header = TRUE, sep = "\t")
row.names(eset_macro) <- eset_macro[, 1]
eset_macro <- eset_macro[, -1]

eset_cd8 <- read.table("/Users/rameshsubramani/Downloads/500_genes/out_CD8_exhausted.tsv", header = TRUE, sep = "\t")
row.names(eset_cd8) <- eset_cd8[, 1]
eset_cd8 <- eset_cd8[, -1]

# Calculate correlation matrices
macro_cor <- cor(t(eset_macro), method="spearman")
cd8_cor <- cor(t(eset_cd8), method="spearman")
dim(macro_cor)
dim(cd8_cor)

# Plotting the correlation matrices
par(mfrow=c(1,2))
plot(im(macro_cor[nrow(macro_cor):1,]),col=cm.colors(256), main="Macro")
plot(im(cd8_cor[nrow(cd8_cor):1,]),col=cm.colors(256), main="CD_8")
#################################################################################################

# Calculate Co-expression matrix for each dataset
cor_matrix_macro <- cor(macro_cor)
macro_graph <- graph.adjacency(cor_matrix_macro, weighted=TRUE, mode="lower")
macro_graph <- delete.edges(macro_graph, E(macro_graph)[ weight < 0.95 ])
macro_graph <- igraph::simplify(macro_graph, remove.multiple = TRUE, remove.loops = TRUE)
write.graph(macro_graph, "macro1_graph.gml", format="gml")

cor_matrix_cd8 <- cor(cd8_cor)
cd8_graph <- graph.adjacency(cor_matrix_cd8, weighted=TRUE, mode="lower")
cd8_graph <- delete.edges(cd8_graph, E(cd8_graph)[ weight < 0.95 ])
cd8_graph <- igraph::simplify(cd8_graph, remove.multiple = TRUE, remove.loops = TRUE)
write.graph(cd8_graph, "cd81_graph.gml", format="gml")
g1 <- macro_graph
g2 <- cd8_graph

# Visualize networks
plot(macro_graph, main = "Co-expression Network for Macrophages")
plot(cd8_graph, main = "Co-expression Network for CD8 Exhausted")

g1.fc <- fastgreedy.community(g1)
sizes(g1.fc)

g2.fc <- fastgreedy.community(g2)
sizes(g2.fc)

hc.macro <- cluster.molecule(eset_macro,method="pearson", linkage="average")  
hc.cd8 <- cluster.molecule(eset_cd8,method="pearson", linkage="average")
gg.macro <- cutree(hc.macro, h=0.4) 
gg.cd8 <- cutree(hc.cd8, h=0.4)
table(gg.macro[table(gg.macro)!=1])
table(gg.macro)
table(gg.cd8)
res.macro <- get.eigen.molecule(eset_macro, groups=gg.macro,  methods="svd", n=2)
res.cd8 <- get.eigen.molecule(eset_cd8, groups=gg.cd8, methods="svd", n=2)


macro_g <- get.eigen.molecule.graph(res.macro)
plot(macro_g, layout=layout.fruchterman.reingold(macro_g), main = "Macro")
write.modules(gg.macro, res.macro, outfile="/Users/rameshsubramani/Downloads/DiffCorr/GSE/module1_list.txt")

cd8_g <- get.eigen.molecule.graph(res.cd8)
plot(cd8_g, layout=layout.fruchterman.reingold(cd8_g), main="CD8") 
write.modules(gg.cd8, res.cd8, outfile="/Users/rameshsubramani/Downloads/DiffCorr/GSE/module2_list.txt")

data2 <- cbind(eset_macro, eset_cd8)
dim(data)

plotDiffCorrGroup(data, gg.macro, gg.cd8, 1, 4, 1:5, 6:7, scale.center=TRUE, scale.scale=TRUE,ylim=c(-5,5))
??plotDiffCorrGroup
comp.2.cc.fdr(output.file="/Users/rameshsubramani/Downloads/DiffCorr/GSE/DiffCorr_res.txt",data[,1:100], data[,101:200], threshold=0.05)
