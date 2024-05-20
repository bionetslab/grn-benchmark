library(optparse)
library(data.table)
library(tidyverse)

# load dcloc and dcglob function
source('12918_2013_1140_MOESM1_ESM.r')

# parsing arguments --------------
option_list <- list(
  make_option(c("-c", "--input.file.1"), type = 'character',
              help = 'Input file'),
  make_option(c("-d", "--input.file.2"), type = 'character',
              help = 'Input file'),
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
  ),
  make_option(
    "--local",
    action = "store_true",
    default = FALSE,
    help = "Run local topology analysis, default is FALSE"
  ),
  make_option(
    "--global",
    action = "store_true",
    default = TRUE,
    help = "Run global topology analysis, default is TRUE"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

print("Starting generation of differential GRN ...")




# functions ------------------------


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
    network.filtered <- cor[df$keepGene, ]
    network.filtered <- network.filtered[, df$keepGene]
    
    
    #### filter correlations over threshold
    
    mask <- (abs(network.filtered) >= corr.threshold)
    network.filtered[!mask] <- 0
    
    return(list(correlation = network.filtered, names = gene.names))
    
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
    network.filtered <- cor[df$keepGene, ]
    network.filtered <- network.filtered[, df$keepGene]
    
    
    #### filter correlations over threshold
    
    mask <- (abs(network.filtered) >= corr.threshold)
    network.filtered[!mask] <- 0
    
    return(list(correlation = network.filtered, names = gene.names))
    
  }




# append given correlation network to output network
merge.network <- function(network, df.corr, condition) {
  # Iterate over each pair of nodes
  for (i in 1:nrow(df.corr)) {
    for (j in 1:ncol(df.corr)) {
      # Skip self-correlations
      if (i != j) {
        # Extract correlation value
        correlation <- df.corr[i, j]
        
        if (correlation != 0) {
          # if correlation exists btw nodes add to network
          # Store the edge in the dataframe
          network <-
            rbind(
              network,
              data.frame(
                target = colnames(df.corr)[i],
                regulator = colnames(df.corr)[j],
                condition = condition,
                weight = correlation,
                edge.type = "undirected"
              )
            )
          
        }
      }
    }
  }
  
  return(network)
}


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


#### local topology analysis

if (opt$local == TRUE) {
  print(paste("Selected topological dissimilarity threshold:", opt$dthresh))
  print(paste("Selected correlation threshold", opt$corrthresh))
  
  print('Starting local topology analysis ...')
  
  df.dcloc = DCloc(data.1[, -1], data.2[, -1])
  df.dcloc <- as.data.frame(df.dcloc)
  
  
  ##### add gene names
  
  df.dcloc$Gene <- data.2$Gene
  df.dcloc <- df.dcloc %>% relocate(Gene)
  
  print("Saving output of DCloc ....")
  
  write.table(
    df.dcloc,
    file = paste(paste0(opt$output.path, "/dcloc_output.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  network.output <-
    network.local(
      df.dcloc,
      1,
      cor.1,
      d.threshold = opt$dthresh,
      corr.threshold = opt$corrthresh
    )
  cond.1.corr <- as.data.frame(network.output$correlation)
  colnames(cond.1.corr) <- network.output$names
  rownames(cond.1.corr) <- network.output$names
  
  network.output.2 <-
    network.local(
      df.dcloc,
      2,
      cor.2,
      d.threshold = opt$dthresh,
      corr.threshold = opt$corrthresh
    )
  cond.2.corr <- as.data.frame(network.output.2$correlation)
  colnames(cond.2.corr) <- network.output.2$names
  rownames(cond.2.corr) <- network.output.2$names
  
  
  
  #### global topology analysis
} else {
  print(paste("Selected p-value threshold:", opt$pthresh))
  print(paste("Selected correlation threshold", opt$corrthresh))
  
  print('Starting global topology analysis ...')
  
  df.dcglob = DCglob(data.1[, -1], data.2[, -1])
  df.dcglob <- as.data.frame(df.dcglob)
  df.dcglob$Gene <- data.1$Gene
  df.dcglob <- df.dcglob %>% relocate(Gene)
  
  
  print("Saving output of DCglob ....")
  
  write.table(
    df.dcglob,
    file = paste(paste0(opt$output.path, "/dcglob_output.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  network.output <-
    network.global(
      df.dcglob,
      1,
      cor.1,
      p.threshold = opt$pthresh,
      corr.threshold = opt$corrthresh
    )
  cond.1.corr <- as.data.frame(network.output$correlation)
  colnames(cond.1.corr) <- network.output$names
  rownames(cond.1.corr) <- network.output$names
  
  network.output.2 <-
    network.global(
      df.dcglob,
      2,
      cor.2,
      p.threshold = opt$pthresh,
      corr.threshold = opt$corrthresh
    )
  cond.2.corr <- as.data.frame(network.output.2$correlation)
  colnames(cond.2.corr) <- network.output.2$names
  rownames(cond.2.corr) <- network.output.2$names
  
}

print("Creating network ...")

# Create an empty dataframe to store diffGRN
network <- data.frame(
  target = character(),
  regulator = character(),
  condition = character(),
  weight = numeric(),
  edge.type = character(),
  stringsAsFactors = FALSE
)


network <- merge.network(network, cond.1.corr, cond.name.1)


network <- merge.network(network, cond.2.corr, cond.name.2)



print("Saving network ...")
fwrite(network, paste0(opt$output.path, "/network.txt"), sep = "\t")
print("Done!")
print(paste("Network saved under", paste0(opt$output.path, "/network.txt")))
