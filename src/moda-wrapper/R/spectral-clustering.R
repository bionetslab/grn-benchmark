#Community Detection using MODA

#Parameter defaults
indicator <-  'X'     # indicator for gene data profile 1
k_clusters <- 5 #Default number of clusters

#Get input args
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<2){
  cat ("few command arguments specified")
  q("no")
}


#extract script input arguments
input_path <- args[1]
output_path <- args[2]
#extract parameters
indicator_arg <- "-i1"
k_clusters_arg <- "-k"

if(indicator_arg %in% args){
  indicator_index <- which(args == indicator_arg) + 1
  indicator <-  args[indicator_index] 
}
# k clusters
if(k_clusters_arg %in% args){
  k_index <- which(args == k_clusters_arg) + 1
  k_clusters <- as.numeric(args[k_index])
}


require(tidyverse)
require(data.table)
require(MODA)


#Print Parameters for this session
print("Setup:")

print(paste0("Label for Gene expression data: ", indicator))
print(paste0("k clusters: ", k_clusters))



output_path_vec <- c('output-spectral') #output directory

for (path in output_path_vec) {
  if(!dir.exists(path)){
    dir.create(paste0(output_path, path))
  }
}



#Load transform
gene_expression_data_path_X=list.files(input_path)[1];



gene_expression_data_X <-fread(file = paste0(input_path, gene_expression_data_path_X))


#transpose data
t_gene_expression_data_X <- t(gene_expression_data_X)
#make a numeric matrix
gene_expr_matrix_X <- as.matrix(t_gene_expression_data_X)

#Spectral Clustering
print("Performing Spectral Clustering...")
spectral_clustering_path <- output_path_vec[1]
ResultFolder = paste0(output_path, spectral_clustering_path)  # where middle files are stored

   
GeneNames <- gene_expression_data_X$Gene
WeightedModulePartitionSpectral(gene_expr_matrix_X[-1,-1],ResultFolder,indicator,
                                GeneNames,k=5) # k is the number of clusters


