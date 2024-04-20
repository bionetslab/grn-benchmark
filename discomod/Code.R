install.packages("devtools")
devtools::install_github("arbet003/discoMod")
install.packages("mclust")
install.packages("foreach")
install.packages("BiocManager")
BiocManager::install("multtest")

library(foreach)
library(devtools)
library(multtest)
library(mclust)

data(golub)

data(golub.cl)
golub.cl

X1 = golub[,which(golub.cl==0)] #ALL
X2 = golub[,which(golub.cl==1)] #AML
rownames(X1) = golub.gnames[,3]
rownames(X2) = golub.gnames[,3]

set.seed(1234)
ind = sample(1:nrow(X1),200)
X1 = X1[ind,]
X2 = X2[ind,]

modules = find_modules(X1,X2,cluster_group=1)
modules$num_modules

ngm = unlist(lapply(modules$group1_modules, ncol))
summary(ngm)

testmods = test_modules(group1_modules = modules$group1_modules, group2_modules = modules$group2_modules)
#View(testmods$pvalues)
#View(testmods$qvalues)

which(testmods$pvalues$PND6 <= 0.05)
which(testmods$qvalues$PND6 <= 0.05)


heat = corrheatmap(modules$group1_modules[[5]], modules$group2_modules[[5]])
plot(heat$both_plots)
plot(heat$group1_plot)
plot(heat$group2_plot)

