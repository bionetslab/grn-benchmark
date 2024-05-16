# Enhancing Gene Network Analysis: A Comparative Study of DiNA Using the p-Norm Difference Test

This study evaluates different statistical methods to identify Differentially Co-expressed Modules (DCMs) in gene networks using Differential Network Analysis (DiNA). Then authors tested a new method, the p-norm difference test (PND), against traditional approaches through simulations and real data analysis. The outcome is the `discoMod` R package, equipped with tools for visualizing networks, testing for DCMs, and clustering genes. Our findings show that PND is an effective alternative to standard methods, as it successfully lowers the false positive rate while keeping a good true positive rate.


 ## Methodology

The study describes a structured set of techniques for performing Differential Network Analysis (DiNA) to find differentially co-expressed modules (DCMs) in gene expression data.


The identification of gene modules is the first step in the procedure. These modules, each of which needs to have three genes at least, are found using clustering techniques or are taken from databases that already exist. It is noteworthy that modules have the ability to overlap, meaning that a gene may be included in more than one module.


Building a gene expression matrix X(gm) for each module across several groups, such as healthy versus diseased individuals, is the next step once the modules have been identified. These matrices' expression levels can be expressed as continuous values if they are acquired via sequencing technologies, or as integer counts from microarray data.

In order to examine the connections within every module, a similarity matrix S(gm) is constructed for every group. The number of genes in the module determines the dimensions of this symmetric matrix. Depending on the type of data and the particular requirements of the analysis, the elements of the similarity matrix are produced using different measures of correlation, such as partial correlation or mutual information, or Pearson, Spearman, or Kendall.


Comparing the network structures represented by the similarity matrices of various groups is the fundamental component of the process, which is based on hypothesis testing. Through the use of several statistical tests, the study specifically looks to examine if there are significant differences between groups in the interactions inside a module between matrices S(g1m) and S(g2m).

The p-norm difference test (PND), a novel strategy the researchers presented, is highlighted as one of the statistical tests used in the publication. With a low false positive rate, this test is intended to efficiently detect notable variations in network topologies.

A non-parametric permutation method is used to compute p-values, which is an essential part of statistical analysis. Because it considers the intricate connections present in the data, this approach is especially well-suited to the structured similarity matrices of the data.


Ultimately, a R package called `discoMod` has the implementation of the full technique along with the tests that go along with it. This program facilitates a streamlined method for researchers by offering extensive capabilities for testing for differential co-expression, grouping gene modules, and displaying the results.

### Hyperparameters:
The parameters and settings for the simulations include:

**Number of Genes in a Module (N)**: The simulations considered modules of different sizes, specifically modules with 10, 50, or 100 genes, depending on the scenario being tested​​.  
**Proportion of Changed Correlations (γ)**: The simulations involved changing a specific proportion of the correlations within the modules to explore different degrees of network alteration. This proportion, denoted as γ, was set at different levels (0.1, 0.4, or 0.7) to represent small, medium, and large effects. In different scenarios, the changes included dropping correlations to zero or altering them by increasing or decreasing their values​​.  
**The compound symmetric correlation parameter (𝜌)**: quantifies the strength of the correlation between any two genes within the module.  
**The p-value of the p-norm difference test (PND)**: calculated using a non-parametric permutation approach, which assesses the significance of the observed differences in network structures between groups by comparing them against a distribution of differences generated under the null hypothesis.

The parameters and settings for the testing serve as foundational structures for analyzing and visualizing the interactions between genes. such as:

**Correlation Matrix**: Measures the similarity between gene expression patterns, with each element representing the correlation coefficient between gene pairs.   
**Adjacency Matrix**: Converts correlation coefficients into a binary or weighted network structure, indicating the presence and strength of connections between genes.
**Topological Overlap Matrix (TOM)**: Enhances the adjacency matrix by considering not only direct connections between genes but also their shared connections, providing a deeper insight into the network's interconnectedness.

**Default**: Discomod authors suggested he paper uses ρ values of 0.3 or 0.7 in different simulation scenarios. These values represent the strength of correlation within the gene modules in the network structure.  While various values from 4 to 20 are tested, the recommended default value for general use is PND6, though the paper advises exploring other values based on specific dataset characteristics .



## Tutorial & Documentation of Discomod Library

For a tutorial, see the example near the end of the following help page:

``` r
library(discoMod)
?test_modules
```

 ## Installation / Run the Project 
 ### Prerequisites
 - R programming language (Version 3.6 or later recommended)
 - RStudio (For running R scripts)

### Libraries
 Need to install the following R packages:

``` r
install.packages("devtools")
devtools::install_github("arbet003/discoMod")
```

### Additional Libraries
```r
install.packages("mclust")
install.packages("foreach")
install.packages("BiocManager")
BiocManager::install("multtest")
```

and After installation, you can load the package and use its functions as follows: 
```r
library(foreach)
library(devtools)
library(multtest)
library(mclust)
```

### Dataset
For the dataset that the paper used from the BiocManager, need to write the followings:
```r
data(golub.cl)
```

As we are using reference data of CD8 and Microphages, data need to be pre-processed first.  
**Data pre processing**:

1. First row and column deletation
2. Row added to rownames
3. String need to convert to numeric

### Functions

We have 3 functions in the GitHub library such as

`find modules`: among the data, it finds modules (where gene >= 3) and it also considers clustering.  
`Test modules`: for DCM, it tested all modules which found via find modules and gave the specific modules ( which Modules are performed DCM)
now we have all analysis data.  
For the better understanding, we use different visualisation representation such as ggplot, correlation heatmap.

 ## Technology used

- **R**: Used for all statistical computing and graphics.
- ***discoMod R package***: Developed specifically for this project to facilitate the differential network analysis.
- ***GitHub***: Used for version control.
 


## How to Use the Project
To use this project, follow these steps:

1. Install the **discoMod** package as shown above.

2. Prepare your gene expression data in the correct format (see Input File Format Specification below).

3. Run the code that is given in code section using functions, passing data as the argument.

4. Analyze the results using the visualization and summary functions provided by the package.


## Input file format specification

For the data, need to use above dataset section.
The input file should be a TSV file with the following specifications:

- Each row represents a gene .
- Each column represents a sample (except the first column, which contains gene identifiers).
- The first row contains sample identifiers


## Output file format specification
The output from the discoMod package functions will typically be an object containing:

- A list of DCMs identified.
- Statistical significance levels for each module.
- Visualization plots (saved as PNG files).

Camparison of all statistical modules and the new PND module. 

## Implementation of Discomod (Code Execution):

### library Installations:
```r
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
```

### Load Reference Data:
```r
# load the data 
data1 <- read.table("./reference_datasets/500_genes/out_CD8_exhausted.tsv")
data2 <- read.table("./reference_datasets/500_genes/out_Macrophages.tsv")
```
As we are using reference data of CD8 and Microphages, data need to be pre-processed first.  
### Data pre processing:

1. First row and column deletation
2. Row added to rownames

```r
# delete top rows (label names)
data1 <- data1[-1,]
data2 <- data2[-1,]

# separate 1st column and set it as row names
rows = data1[, 1]
data1 = data1[, -1]
rownames(data1) = rows

rows = data2[, 1]
data2 = data2[, -1]
rownames(data2) = rows
```
Among the proceed data, it finds modules (where gene >= 3)
```r
#Find the module by using find_module functions
modules = discoMod::find_modules(data1, data2, cluster_group=2)
modules$num_modules
```

3.  String need to convert to numeric:
```r
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
```
Among all the module, test specific modules ( which Modules are performed DCM):
```r
testmods = discoMod::test_modules(group1_modules = modules$group1_modules , group2_modules = modules$group2_modules)

```
### Visualization using Correlation heatmap & ggplot
```r
heat = discoMod::corrheatmap(modules$group1_modules[[1]], modules$group2_modules[[1]])
#for correlation heatmap
plot(heat$both_plots)
# plot(heat$group1_plot)
# plot(heat$group2_plot)

#for line graph
library(tidyr)
library(viridis)
library(ggplot2)

graph_data = testmods$pvalues
graph_data$module = rownames(graph_data)
graph_data <- pivot_longer(graph_data, cols= -module, names_to = 'Test', values_to = 'Values')

plot <- ggplot(data = graph_data, aes(x = Test, y = log10(Values), group= module, color= module)) +
geom_line() +
geom_point() +
labs(x= 'Test', y= 'Values', title = 'Plot') +
theme_minimal() +
scale_colour_viridis(discrete = TRUE)

plot1 <- plot  + ylim(-10, 3)
plot

```
Reference Dataset:
![full-cluster2-53modules](https://github.com/NaziaUrmi/grn-benchmark/assets/46952253/256d1de4-a170-450a-835f-358177bac1e2)
    
GEO Dataset:

![new-5mod 23genes](https://github.com/NaziaUrmi/grn-benchmark/assets/46952253/36188c8b-7a0d-4608-9fc2-026521019e45)


## Result 

The analysis using the Golub dataset demonstrated that the PND method significantly outperformed traditional statistical approaches in identifying Differentially Co-expressed Modules (DCMs). This was particularly notable in the accuracy of detecting true positive changes in gene network structures, making PND a robust choice for datasets with similar characteristics.  
In contrast, the application of the PND method to the CD8 and Macrophages datasets did not yield as definitive results. The method encountered challenges in consistently identifying significant differences in the network structures. This inconsistency suggests potential limitations of PND when applied to datasets with certain biological complexities or differing underlying network dynamics. 



## Reference

The discoMod R package:
Arbet, J., Zhuang, Y., Litkowski, E., Saba, L., & Kechris, K. (2021). Comparing Statistical Tests for Differential Network Analysis of Gene Modules. Frontiers in genetics, 12, 748.


