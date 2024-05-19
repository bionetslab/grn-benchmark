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
devtools::install_github("arbet003/discoMod")

library(mclust)


data <- read.table("./reference_datasets/GSE206793_expression.tsv")
data1 <- data %>% select(contains("prostate.cancer_M1"))
data2 <- data %>% select(contains("prostate.cancer_N0M0"))
data2 <- data2[, 0:31]

# delete top rows (label names)
#data1 <- data1[-1,]
#data2 <- data2[-1,]

# separate 1st column and set it as row names
#rows = data1[, 1]
#data1 = data1[, -1]
#rownames(data1) = rows

#rows = data2[, 1]
#data2 = data2[, -1]
#rownames(data2) = rows

set.seed(1234)
ind = sample(1:nrow(data1),200)
data1 = data1[ind,]
data2 = data2[ind,]


#Find the module by using find_module functions
modules = discoMod::find_modules(data1, data2, cluster_group=2)
modules$num_modules


# convert to numeric format
for (i in seq_along(modules$group1_modules)) {
  if (length(modules$group1_modules[[i]]) > 0) {
    modules$group1_modules[[i]] = apply(modules$group1_modules[[i]], c(1, 2), as.numeric)
  }
}

for (i in seq_along(modules$group2_modules)) {
  if (length(modules$group2_modules[[i]]) > 0) {
    modules$group2_modules[[i]] = apply(modules$group2_modules[[i]], c(1, 2), as.numeric)
  }
}

# 
# library(parallel)
# library(doSNOW)
# testmods = discoMod::test_modules(group1_modules = modules$group1_modules, 
#                                   group2_modules = modules$group2_modules,
#                                   cortype = "pearson",
#                                   othertype = "Adjacency",
#                                   cores = 8, 
#                                   perm_stats=FALSE)
testmods = discoMod::test_modules(group1_modules = modules$group1_modules , group2_modules = modules$group2_modules)


#getting corelation matrix to proof the arbitary hypothesis

heat = discoMod::corrheatmap(modules$group1_modules[[1]], modules$group2_modules[[1]])
# for heatmap
plot(heat$both_plots)
# plot(heat$group1_plot)
# plot(heat$group2_plot)

#for line graph

library(tidyr)
library(viridis)
library(ggplot2)

graph_data = testmods$qvalues
graph_data$module = rownames(graph_data)
graph_data <- pivot_longer(graph_data, cols= -module, names_to = 'Test', values_to = 'Values')

plot <- ggplot(data = graph_data, aes(x = Test, y = (Values), group= module, color= module)) +
  geom_line() +
  geom_point() +
  labs(x= 'Test', y= 'Values', title = 'Plot') +
  theme_minimal() +
  scale_colour_viridis(discrete = TRUE)


plot
