if (!require('devtools')) {
  install.packages('devtools')
}
if (!require('mclust')) {
  install.packages('mclust')
}
if (!require('foreach')) {
  install.packages('foreach')
}

library(foreach)
library(devtools)
library(multtest)
library(mclust)

#data1 <- read.table("C:/Users/Syed Hussain/bionet/grn-benchmark/reference_datasets/500_genes/out_CD8_exhausted.tsv")
#data2 <- read.table("C:/Users/Syed Hussain/bionet/grn-benchmark/reference_datasets/500_genes/out_Macrophages.tsv")

X1 = golub[,which(golub.cl==0)] #ALL
X2 = golub[,which(golub.cl==1)] #AML
rownames(X1) = golub.gnames[,3]
rownames(X2) = golub.gnames[,3]

set.seed(1234)
ind = sample(1:nrow(X1),200)
X1 = X1[ind,]
X2 = X2[ind,]

# data_matrix1 <- data1[-1, -1]
# data_matrix2 <- data2[-1, -1]
# 
# rownames(data_matrix1) = data1[-1,1]
# rownames(data_matrix2) = data2[-1,1]
# 
# modules = discoMod::find_modules(data_matrix1, data_matrix2, cluster_group=1)
# modules$num_modules

modules = discoMod::find_modules(X1, X2, cluster_group=1)
modules$num_modules

ngm = unlist(lapply(modules$group1_modules, ncol))
summary(ngm)

# Filter out list elements with length zero
filtered_list1 <- modules$group1_modules[sapply(modules$group1_modules, length) > 0]
filtered_list2 <- modules$group2_modules[sapply(modules$group2_modules, length) > 0]

testmods = discoMod::test_modules(group1_modules = modules$group1_modules, group2_modules = modules$group2_modules)

#View(testmods$pvalues)
#View(testmods$qvalues)

which(testmods$pvalues$PND6 <= 0.05)
which(testmods$qvalues$PND6 <= 0.05)


heat = corrheatmap(modules$group1_modules[[5]], modules$group2_modules[[5]])
plot(heat$both_plots)
plot(heat$group1_plot)
plot(heat$group2_plot)

# # Load the ggplot2 package
# library(ggplot2)
# 
# # Load the tidyr library for data manipulation
# library(tidyr)
# 
# graph_data = testmods$test_stats
# graph_data$module <- rownames(graph_data)
# 
# # Reshape the data from wide to long format
# graph_data <- pivot_longer(graph_data, cols = -module, names_to = "Test", values_to = "Values")
# 
# 
# library(viridis)
# 
# 
# # Plotting using ggplot2
# ggplot(data = graph_data, aes(x = Test, y = log10(Values), group = module, color = module)) +
#   geom_line() +  # Add lines to connect points
#   geom_point() +  # Add points
#   labs(x = "Test", y = "Values", title = "Connected Scatter Plot") +
#   theme_minimal() +
#   scale_color_viridis(discrete = TRUE)  # Improves color selection


