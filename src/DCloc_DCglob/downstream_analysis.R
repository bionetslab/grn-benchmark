library(data.table)
library(tidyverse)
library(gplots)
library(optparse)
library(igraph)


# parsing arguments --------------
option_list <- list(
  make_option(c("-c", "--input.file.1"), type = 'character',
              help = 'Input file'),
  make_option(c("-d", "--input.file.2"), type = 'character',
              help = 'Input file'),
  make_option("--dcloc.output", type = 'character',
              help = 'Output file of DCloc'),
  make_option("--dcglob.output", type = 'character',
              help = 'Output file of DCglob'),
  make_option(c("-o", "--output.path"), type = 'character',
              help = 'output path'),
  make_option(
    "--dthresh",
    type = 'numeric',
    default = 0.3,
    help = 'Local topology analysis: threshold for topological dissimilarity'
  ),
  make_option(
    "--pthresh",
    type = 'numeric',
    default = 0.05,
    help = 'Global topology analysis: threshold for p-value'
  ),
  make_option(
    "--corrthresh",
    type = 'numeric',
    default = 0.5,
    help = 'Threshold for correlation matrix for which edges to keep'
  )
)

opt <- parse_args(OptionParser(option_list = option_list))



# functions --------------------------------------------------------------------------------------------------------------------------------------


# creates network topology based on differentially correlated genes of the global topology analysis


# calculate correlation network of differentially correlated genes of global topology analysis


network.global <-
  function(df,
           condition,
           cor,
           p.threshold = 0.1,
           corr.threshold = 0.5) {
    # Parameters:
    # df: dataframe of global topology analysis, dcglob output
    # condition: condition of which the differentially correlated genes should be analyzed
    # p.threshold: threshold to choose cutoff of differentially correlated genes
    # corr.threshold: correlation threshold for correlation network to only include edges with weights higher than corr.thresh
    # cor: correlation matrix of condition
    
    
    # filter genes important for given condition
    df$keepGene.cond <- (df$`High. Corr.` == condition)
    df$keepGene.p <- (df$`p-value` < p.threshold)
    df$keepGene <- df$keepGene.cond & df$keepGene.p
    gene.names <- df$Gene[df$keepGene]
    network.filtered <- cor[df$keepGene,]
    network.filtered <- network.filtered[, df$keepGene]
    
    
    #### filter correlations over threshold
    
    mask <- (abs(network.filtered) >= corr.threshold)
    
    
    return(list(
      correlation = mask,
      names = gene.names,
      df = df
    ))
    
  }

# calculate correlation network of differentially correlated genes of local topology analysis
network.local <-
  function(df,
           condition,
           cor,
           d.threshold = 0.25,
           corr.threshold = 0.5) {
    # Parameters:
    # df: dataframe of local topology analysis; dcloc output
    # condition: condition of which the differentially correlated genes should be analyzed
    # cor: correlation matrix of condition
    # d.threshold: threshold to choose cutoff of differentially correlated genes based on topological dissimilarity
    # corr.threshold: correlation threshold for correlation network to only include edges with weights higher than corr.thresh
    
    
    # filter genes important for given condition
    df$Corr <- ifelse(df$mean.neigh.A > df$mean.neigh.B, 1, 2)
    df$keepGene.cond <- (df$Corr == condition)
    df$keepGene.d <- (df$abs.top.dissim > d.threshold)
    df$keepGene <- df$keepGene.cond & df$keepGene.d
    gene.names <- df$Gene[df$keepGene]
    network.filtered <- cor[df$keepGene,]
    network.filtered <- network.filtered[, df$keepGene]
    
    
    #### filter correlations over threshold
    
    mask <- (abs(network.filtered) >= corr.threshold)
    
    
    return(list(
      correlation = mask,
      names = gene.names,
      df = df
    ))
    
  }


cluster.genes <- function(data, cor.threshold = 0.4) {
  #output cluster.colors: cluster colors for cor.threshold and clustering: output of hierarchical clustering, clusters: clusters with cutoff threshold
  
  
  distance.matrix <- as.dist(1 - cor(t(data)))
  
  # Perform hierarchical clustering
  hc <- hclust(distance.matrix, method = "average")
  
  # Cut the tree at a specified Pearson correlation threshold
  height.threshold <- 1 - cor.threshold
  clusters <- cutree(hc, h = height.threshold)
  
  # Define colors for clusters
  cluster.colors <- rainbow(length(unique(clusters)))[clusters]
  
  #output cluster colors for cor.threshold and output of hierarchical clustering,
  return(list(
    cluster.colors = cluster.colors,
    clustering = hc,
    clusters = clusters
  ))
  
  
}

#saves clustering of given condition in a text file

save_clustering <- function(output.path, condition.name, cluster, row.names, global = TRUE){

  cluster_list <- split(row.names, cluster)
  alg = ifelse(global, "global", "local")
  
  fileConn<-file(paste0(output.path, "/clusters_", alg, "_", condition.name, ".txt"), "a")
  
  # Write clusters to the file
  for (cluster_id in names(cluster_list)) {
    writeLines(paste("Cluster", cluster_id, ":"), fileConn)
    writeLines(cluster_list[[cluster_id]], fileConn)
  }
  
  close(fileConn)
  
}

#function centers data to mean 0 and std 1
center.data <- function(data) {
  centered.data <- t(scale(t(data), center = TRUE, scale = TRUE))
  return(centered.data)
}


#analyzes the expression data of the given conditions based on the stronger correlated genes of given condition to analyze
# creates heatmaps with hierarchical clustering
# creates individual gml graphs that can be loaded in cytoscape

downstream_analysis.local <-
  function(df.dcloc,
           # dataframe with output of dcloc
           condition.to.analyze,
           # either 1 or 2
           cond1.name,
           #name of condition 1
           cond2.name,
           # name of condition 2
           d.threshold,
           # threshold of topological dissimilarity
           corr.threshold,
           #correlation threshold
           cor.1,
           # correlation matrix of condition 1
           cor.2,
           # correlation matrix of condition 1
           data.cond1,
           #expression data of condition 1
           data.cond2,
           # expression data of condition 2
           output.path) {
    # path where to save graphs, heatmaps and clustering list
    
    
    network.1 <-
      network.local(
        df.dcloc,
        condition.to.analyze,
        cor.1,
        d.threshold = d.threshold,
        corr.threshold = corr.threshold
      )
    
    cond.1.corr <- as.data.frame(network.1$correlation)
    
    colnames(cond.1.corr) <- network.1$names
    rownames(cond.1.corr) <- network.1$names
    
    network.2 <-
      network.local(
        df.dcloc,
        condition.to.analyze,
        cor.2,
        d.threshold = d.threshold,
        corr.threshold = corr.threshold
      )
    
    cond.2.corr <- as.data.frame(network.2$correlation)
    colnames(cond.2.corr) <- network.2$names
    rownames(cond.2.corr) <- network.2$names
    
    graph.1 <-
      graph.adjacency(network.1$correlation, mode = c("undirected"))
    graph.2 <-
      graph.adjacency(network.2$correlation, mode = c("undirected"))
    
    gene.names <- network.2$names
    
    #### assign gene names
    V(graph.1)$name <- gene.names
    V(graph.2)$name <- gene.names
    
    cond.name = ifelse(condition.to.analyze == 1, cond1.name, cond2.name)
    
    ### save graph
    save.1 <-
      paste0(output.path,
             "/graph_DCloc_",
             cond.name,
             "_",
             cond1.name,
             ".gml")
    save.2 <-
      paste0(output.path,
             "/graph_DCloc_",
             cond.name,
             "_",
             cond2.name,
             ".gml")
    
    write.graph(graph.1, save.1, format = "gml")
    write.graph(graph.2, save.2, format = "gml")
    
    
    
    ###### plot heatmap
    
    
    #filter out relevant genes
    rel.data.cond1 <- data.cond1[network.1$df$keepGene, -1]
    rel.data.cond2 <- data.cond2[network.1$df$keepGene, -1]
    
    rownames(rel.data.cond1) <- gene.names
    rownames(rel.data.cond2) <- gene.names
    
    
    data.centered.data.cond2 <- center.data(rel.data.cond2)
    data.centered.data.cond1 <- center.data(rel.data.cond1)
    
    cluster.1 <- cluster.genes(data.centered.data.cond1)
    cluster.2 <- cluster.genes(data.centered.data.cond2)
    
    
    if(condition.to.analyze == 1){
      cluster <- cluster.1
    } else {
      cluster <- cluster.2
    }
    
    
    save_clustering(
      output.path = output.path,
      cond.name,
      cluster.1$clusters,
      rownames(rel.data.cond1),
      global = FALSE
    )
    
    
    pdf(paste0(output.path, "/Heatmap_DCloc", cond.name, "_", cond1.name, ".pdf"))
    
    heatmap(
      as.matrix(rel.data.cond1),
      RowSideColors = cluster$cluster.colors,
      Rowv = as.dendrogram(cluster.1$clustering),
      labRow = rownames(rel.data.cond1),
      labCol = colnames(rel.data.cond1),
      main = cond1.name,
      cexRow = 0.2,
      cexCol = 0.2,
      margins = c(9, 3),
      xlab = "Samples",
      ylab = paste0("Differentially correlated genes of ", cond.name)
    )
    
    dev.off()
    
    pdf(paste0(output.path, "/Heatmap_DCloc", cond.name, "_", cond2.name, ".pdf"))
    
    heatmap(
      as.matrix(rel.data.cond2),
      RowSideColors = cluster$cluster.colors,
      Rowv = as.dendrogram(cluster.2$clustering),
      labRow = rownames(rel.data.cond1),
      labCol = colnames(rel.data.cond1),
      main = cond2.name,
      cexRow = 0.2,
      cexCol = 0.2,
      margins = c(9, 3),
      xlab = "Samples",
      ylab = paste0("Differentially correlated genes of ", cond.name)
    )
    
    dev.off()
    
    
    
    
  }



# analyzes the expression data of the given conditions based on the stronger correlated genes of given condition to analyze
# creates heatmaps with hierarchical clustering
# creates individual gml graphs that can be loaded in cytoscape

downstream_analysis.global <-
  function(df.dcglob,
           # dataframe with output of dcloc
           condition.to.analyze,
           # either 1 or 2
           cond1.name,
           #name of condition 1
           cond2.name,
           # name of condition 2
           p.threshold,
           # threshold of topological dissimilarity
           corr.threshold,
           #correlation threshold
           cor.1,
           # correlation matrix of condition 1
           cor.2,
           # correlation matrix of condition 1
           data.cond1,
           #expression data of condition 1
           data.cond2,
           # expression data of condition 2
           output.path) {
    # path where to save graphs, heatmaps and clustering list
    
    
    
    network.1 <-
      network.global(
        df.dcglob,
        condition.to.analyze,
        cor.1,
        p.threshold = p.threshold,
        corr.threshold = corr.threshold
      )
    
    
    cond.1.corr <- as.data.frame(network.1$correlation)
    colnames(cond.1.corr) <- network.1$names
    rownames(cond.1.corr) <- network.1$names
    
    network.2 <-
      network.global(
        df.dcglob,
        condition.to.analyze,
        cor.2,
        p.threshold = p.threshold,
        corr.threshold = corr.threshold
      )
    
    cond.2.corr <- as.data.frame(network.2$correlation)
    colnames(cond.2.corr) <- network.2$names
    rownames(cond.2.corr) <- network.2$names
    
    graph.1 <-
      graph.adjacency(network.1$correlation, mode = c("undirected"))
    graph.2 <-
      graph.adjacency(network.2$correlation, mode = c("undirected"))
    
    gene.names <- network.2$names
   
    
    
    #### assign gene names
    V(graph.1)$name <- gene.names
    V(graph.2)$name <- gene.names
   
    cond.name = ifelse(condition.to.analyze == 1, cond1.name, cond2.name)
    
    ### save graph
    save.1 <-
      paste0(output.path,
             "/graph_DCglob_",
             cond.name,
             "_",
             cond1.name,
             ".gml")
    save.2 <-
      paste0(output.path,
             "/graph_DCglob_",
             cond.name,
             "_",
             cond2.name,
             ".gml")
    
    write.graph(graph.1, save.1, format = "gml")
    write.graph(graph.2, save.2, format = "gml")
    
    
    
    ###### plot heatmap
    
    
    #filter out relevant genes
    rel.data.cond1 <- data.cond1[network.1$df$keepGene, -1]
    rel.data.cond2 <- data.cond2[network.1$df$keepGene, -1]
    

    
    rownames(rel.data.cond1) <- gene.names
    rownames(rel.data.cond2) <- gene.names
    
    data.centered.data.cond2 <- center.data(rel.data.cond2)
    data.centered.data.cond1 <- center.data(rel.data.cond1)
    
    cluster.1 <- cluster.genes(data.centered.data.cond1)
    cluster.2 <- cluster.genes(data.centered.data.cond2)

    if(condition.to.analyze == 1){
      cluster <- cluster.1
    } else {
      cluster <- cluster.2
    }
    


    save_clustering(
      output.path = output.path,
      cond.name,
      cluster$clusters,
      rownames(rel.data.cond1),
      global = TRUE
    )
    
   
    
    pdf(paste0(output.path, "/Heatmap_DCglob", cond.name, "_", cond1.name, ".pdf"))
    
    heatmap(
      as.matrix(rel.data.cond1),
      RowSideColors = cluster$cluster.colors,
      Rowv = as.dendrogram(cluster.1$clustering),
      labRow = rownames(rel.data.cond1),
      labCol = colnames(rel.data.cond1),
      main = cond1.name,
      cexRow = 0.2,
      cexCol = 0.2,
      margins = c(9, 3),
      xlab = "Samples",
      ylab = paste0("Differentially correlated genes of ", cond.name)
    )
    
    
    dev.off()
    
    
    
    pdf(paste0(output.path, "/Heatmap_DCglob", cond.name, "_", cond2.name, ".pdf"))
    
    heatmap(
      as.matrix(rel.data.cond2),
      RowSideColors = cluster$cluster.colors,
      Rowv = as.dendrogram(cluster.2$clustering),
      labRow = rownames(rel.data.cond1),
      labCol = colnames(rel.data.cond1),
      main = cond2.name,
      cexRow = 0.2,
      cexCol = 0.2,
      margins = c(9, 3),
      xlab = "Samples",
      ylab = paste0("Differentially correlated genes of ", cond.name)
    )
    
    dev.off()
    
    
    
    
    
  }




# getting condition names --------
# copied from zscore


data.1 <- fread(file = opt$input.file.1,
                header = TRUE,
                sep = "\t")
data.2 <- fread(file = opt$input.file.2,
                header = TRUE,
                sep = "\t")

# getting condition names --------
# copied from zscore

cond.name.1 <- sub('.*out_', "", opt$input.file.1)
cond.name.2 <- sub('.*out_', "", opt$input.file.2)

cond.name.1 <- sub(".tsv", "", cond.name.1)
cond.name.2 <- sub(".tsv", "", cond.name.2)

data.1 <- fread(file = opt$input.file.1,
                header = TRUE,
                sep = "\t")
data.2 <- fread(file = opt$input.file.2,
                header = TRUE,
                sep = "\t")



##### calculate correlation matrix
Mat.1 <- data.1[, -1]
Mat.2 <- data.2[, -1]
ngenes <- nrow(Mat.1)
cor.1 <- cor(t(Mat.1))
cor.2 <- cor(t(Mat.2))
for (i in 1:ngenes) {
  cor.1[i, i] <- 0
  cor.2[i, i] <- 0
}


dcglob <-  fread(file = opt$dcglob.output, sep = "\t")
dcloc <- fread(file = opt$dcloc.output, sep = "\t")


##### convert to data frame
df.dcglob <- as.data.frame(dcglob)
df.dcloc <- as.data.frame(dcloc)

##### add gene names
df.dcglob$Gene <- data.1$Gene
df.dcloc$Gene <- data.1$Gene
df.dcglob <- df.dcglob %>% relocate(Gene)
df.dcloc <- df.dcloc %>% relocate(Gene)



downstream_analysis.global(
  df.dcglob = df.dcglob,
  condition.to.analyze = 1,
  cond1.name = cond.name.1,
  cond2.name = cond.name.2,
  p.threshold = opt$pthresh,
  corr.threshold = opt$corrthresh,
  cor.1 = cor.1,
  cor.2 = cor.2,
  data.cond1 = data.1,
  data.cond2 = data.2,
  output.path = opt$output.path
)

downstream_analysis.global(
  df.dcglob = df.dcglob,
  condition.to.analyze = 2,
  cond1.name = cond.name.1,
  cond2.name = cond.name.2,
  p.threshold = opt$pthresh,
  corr.threshold = opt$corrthresh,
  cor.1 = cor.1,
  cor.2 = cor.2,
  data.cond1 = data.1,
  data.cond2 = data.2,
  output.path = opt$output.path
)


downstream_analysis.local(
  df.dcloc = df.dcloc,
  condition.to.analyze = 1,
  cond1.name = cond.name.1,
  cond2.name = cond.name.2,
  d.threshold = opt$dthresh,
  corr.threshold = opt$corrthresh,
  cor.1 = cor.1,
  cor.2 = cor.2,
  data.cond1 = data.1,
  data.cond2 = data.2,
  output.path = opt$output.path
)
downstream_analysis.local(
  df.dcloc = df.dcloc,
  condition.to.analyze = 2,
  cond1.name = cond.name.1,
  cond2.name = cond.name.2,
  d.threshold = opt$dthresh,
  corr.threshold = opt$corrthresh,
  cor.1 = cor.1,
  cor.2 = cor.2,
  data.cond1 = data.1,
  data.cond2 = data.2,
  output.path = opt$output.path
)

