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

data1 <- cbind(eset_macro, eset_macro) 
dim(data1)

hc.macro <- cluster.molecule(data1[1:100],method="pearson", linkage="average")  
hc.cd8 <- cluster.molecule(data1[101:200],method="pearson", linkage="average")
gg.macro <- cutree(hc.macro, h=0.4) 
gg.cd8 <- cutree(hc.cd8, h=0.4)
table(gg.macro[table(gg.macro)!=1])
table(gg.macro)
table(gg.cd8)
res.macro <- get.eigen.molecule(data1, groups=gg.macro, methods="svd", n=2)
res.cd8 <- get.eigen.molecule(data1, groups=gg.cd8, methods="svd", n=2)


macro_g <- get.eigen.molecule.graph(res.macro)
plot(macro_g, layout=layout.fruchterman.reingold(macro_g), main = "Macro")
write.modules(gg.macro, res.macro, outfile="module1_list.txt")

cd8_g <- get.eigen.molecule.graph(res.cd8)
plot(cd8_g, layout=layout.fruchterman.reingold(cd8_g), main="CD8") 
write.modules(gg.cd8, res.cd8, outfile="module2_list.txt")

data2 <- cbind(eset_macro, eset_cd8)

plotDiffCorrGroup(data2, gg.macro, gg.cd8, 1, 4, 1:5, 6:7, scale.center=TRUE, scale.scale=TRUE,ylim=c(-5,5))
comp.2.cc.fdr(output.file="/Users/rameshsubramani/Downloads/DiffCorr/GSE/DiffCorr_res.txt",data2[,1:100], data2[,101:200], threshold=0.05)
