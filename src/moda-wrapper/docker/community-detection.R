#Community Detection using MODA

#Parameter defaults
indicator1 <-  'X'     # indicator for gene data profile 1

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
indicator1_arg <- "-i1"

if(indicator1_arg %in% args){
  indicator1_index <- which(args == indicator1_arg) + 1
  indicator1 <-  args[indicator1_index] 
}


require(tidyverse)
require(data.table)
require(MODA)

#Print Parameters for this session
print("Setup:")

print(paste0("Label for Gene expression data: ", indicator1))

output_path_vec <- c('output-louvain') #output directory

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

#Community Detection
print("Performing Community Detection...")
community_detection_path <- output_path_vec[1]
ResultFolder = paste0(output_path, community_detection_path)  # where middle files are stored

GeneNames <- gene_expression_data_X$Gene
intModules <- WeightedModulePartitionLouvain(gene_expr_matrix_X[-1, -1],ResultFolder,indicator1,GeneNames)



