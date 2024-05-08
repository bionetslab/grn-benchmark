# Discomod: Comparing Statistical Tests for Differential Network Analysis of Gene Modules

This study evaluates different statistical methods to identify Differentially Co-expressed Modules (DCMs) in gene networks using Differential Network Analysis (DiNA). Then authors tested a new method, the p-norm difference test (PND), against traditional approaches through simulations and real data analysis. The outcome is the `discoMod` R package, equipped with tools for visualizing networks, testing for DCMs, and clustering genes. Our findings show that PND is an effective alternative to standard methods, as it successfully lowers the false positive rate while keeping a good true positive rate.

 ## Documentation
 This repository contains the implementation of the study on Differential Network Analysis (DiNA) to identify Differentially Co-expressed Modules (DCMs) using various statistical tests, including the newly proposed p-norm difference test (PND). The study involves extensive simulations and analyses on real-world gene expression data to evaluate and compare the effectiveness of different methods.

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
Data pre processing:
- first row and column deletation 
- Row added to rownames
- String need to convert to numeric 

We have 3 functions in the GitHub library such as 
1. find modules: among the data, it finds modules (where gene >= 3) and it also considers clustering 
2. Test modules: for DCM, it tested all modules which found via find modules and gave the specific modules ( which Modules are performed DCM)
3. now we have all analysis data. For the better understanding, we use different visualisation representation such as ggplot, correlation heatmap. 


 ## technology used

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

## Result 

The p-norm difference test (PND) introduced in this study showed superior performance in identifying DCMs with a high true positive rate and a controlled false positive rate compared to existing methods. These results suggest that PND is a robust tool for analyzing complex gene expression data, particularly useful in studies involving diseases or developmental stages.

## Tutorial

For a tutorial, see the example near the end of the following help page:

``` r
library(discoMod)
?test_modules
```

## Reference

When citing the discoMod R package, please use the following:

Arbet, J., Zhuang, Y., Litkowski, E., Saba, L., & Kechris, K. (2021). Comparing Statistical Tests for Differential Network Analysis of Gene Modules. Frontiers in genetics, 12, 748.


