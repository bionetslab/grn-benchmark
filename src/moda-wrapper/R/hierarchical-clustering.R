#Hierarchical Clustering using MODA
#Const
vec_cutting_criterion <- c('Density', 'Modularity')

#Parameter defaults
indicator1 <-  'X'     # indicator for gene data profile 1
indicator2 <-  'Y'      # indicator for gene data profile 2
CuttingCriterion = 'Density' # could be Density or Modularity
specificTheta = 0.1 #threshold to define condition specific modules
conservedTheta = 0.1#threshold to define conserved modules

#Get input args
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<3){
  cat ("few command arguments specified")
  q("no")
}

#extract script input arguments
input_path_1 <- args[1]
input_path_2 <- args[2]
output_path <- args[3]
#extract parameters
indicator1_arg <- "-i1"
indicator2_arg <- "-i2"
cutting_criterion_arg <- "-c"
specific_theta_arg <- "-s"
conserved_theta_arg <- "-t"

if(indicator1_arg %in% args){
  indicator1_index <- which(args == indicator1_arg) + 1
  indicator1 <-  args[indicator1_index] 
}

if(indicator2_arg %in% args){
  indicator2_index <- which(args == indicator2_arg) + 1
  indicator2 <-  args[indicator2_index] 
  
}
if(cutting_criterion_arg %in% args){
  cutting_criterion_index <- which(args == cutting_criterion_arg) + 1
  c_value <-  as.numeric(args[cutting_criterion_index])
  #print(paste0("C_VALUE: ", class(c_value)))
  CuttingCriterion <-  vec_cutting_criterion[c_value] # [1,2] -> [Density, Modularity]
}

if(specific_theta_arg %in% args){
  specific_theta_index <- which(args == specific_theta_arg) + 1
  specificTheta <- as.numeric(args[specific_theta_index])
}
if(conserved_theta_arg %in% args){
  conserved_theta_index <- which(args == conserved_theta_arg) + 1
  conservedTheta <- as.numeric(args[conserved_theta_index])
}



require(tidyverse)
require(data.table)
require(MODA)


#Print Parameters for this session
print("Setup:")

print(paste0("Label for Condition 1: ", indicator1))
print(paste0("Labelfor Condition 2: ", indicator2))
print(paste0("Cutting Criterion: ", CuttingCriterion))
print(paste0("Specific  Theta: ", specificTheta))
print(paste0("Specific  Theta: ", conservedTheta))

output_path_vec <- c('output-hierarchical') #output directory

for (path in output_path_vec) {
  if(!dir.exists(path)){
    dir.create(paste0(output_path, path))
  }
}



#Load transform
gene_expression_data_path_X=list.files(input_path_1)[1];
gene_expression_data_path_Y=list.files(input_path_2)[1];


gene_expression_data_X <-fread(file = paste0(input_path_1, gene_expression_data_path_X))
gene_expression_data_Y <-fread(file = paste0(input_path_2, gene_expression_data_path_Y))

#transpose data
t_gene_expression_data_X <- t(gene_expression_data_X)
t_gene_expression_data_Y <- t(gene_expression_data_Y)
#make a numeric matrix
gene_expr_matrix_X <- as.matrix(t_gene_expression_data_X)
gene_expr_matrix_Y <- as.matrix(t_gene_expression_data_Y)

print("Performing Hierachical Clustering...")
hierachical_clustering_path <- output_path_vec[1]
ResultFolder = paste0(output_path, hierachical_clustering_path, "/")  # where middle files are stored

##modules detection for network 1
intModules1 <- WeightedModulePartitionHierarchical(gene_expr_matrix_X[-1,-1],ResultFolder,
                                                   indicator1,CuttingCriterion)
##modules detection for network 2
intModules2 <- WeightedModulePartitionHierarchical(gene_expr_matrix_Y[-1, -1],ResultFolder,
                                                   indicator2,CuttingCriterion)

#heatmap of correlation matrix of gene expression profile 1
png(paste0(ResultFolder, paste0('/heatdata_', indicator1, '.png')))
heatmap(cor(as.matrix(gene_expression_data_X[-1,-1])))
#dev.off()

#heatmap of correlation matrix of gene expression profile 2
png(paste0(ResultFolder, paste0('/heatdata_', indicator2, '.png')))
heatmap(cor(as.matrix(gene_expression_data_Y[-1,-1])))
dev.off()


#Perform Fix and write identified module gene names
DenseModuleGeneID_files <- list.files(ResultFolder, pattern = "^DenseModuleGeneID_")
DenseModuleGene_files <- character(length(DenseModuleGeneID_files))

counter <- 0
print("Fixing Gene Names in Modules")
cat("\n\n")
for (f in DenseModuleGeneID_files) {
  counter <- counter + 1
  filename_tail <- gsub("DenseModuleGeneID_","", f)
  DenseModuleGene_files[counter] <- paste0("DenseModuleGene_", filename_tail)
  #extract TCells in TCells_1.txt
  gene_expression_indicator <- gsub("^(.*?)_\\d+\\.txt$", "\\1", filename_tail) # extract the indicator name. eg: X or Y
  if(gene_expression_indicator==indicator1){
    module_genesIDs <- fread(file = paste0(ResultFolder, DenseModuleGeneID_files[counter]))
    #ids <- scan(paste0(ResultFolder, DenseModuleGeneID_files[counter]), what = numeric(), sep = "\n")
    #Load Gene from module file
    gene_names_vec=c()
    for (id in module_genesIDs) {
      gene_names_vec <- append(gene_names_vec, gene_expr_matrix_X[1, id])
    }
    
    #write gene names back to file
    writeLines(gene_names_vec, paste0(ResultFolder, DenseModuleGene_files[counter]))
    
  }else{
    #indicator2
    module_genesIDs <- fread(file = paste0(ResultFolder, DenseModuleGeneID_files[counter]))
    #Load Gene from module file
    gene_names_vec=c()
    for (id in module_genesIDs) {
      gene_names_vec <- append(gene_names_vec, gene_expr_matrix_Y[1, id])
    }
    #write gene names back to file
    writeLines(gene_names_vec, paste0(ResultFolder, DenseModuleGene_files[counter]))
  }
  
  
  print(paste0(f," : ", DenseModuleGene_files[counter]))

}



