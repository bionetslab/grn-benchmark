# Discordant
***Version: 1.28.0***

## Input:
- -omics data across two conditions (e.g., gene expression data across two conditions) as a **Bioconductor's ExpressionSet** object.
  - one -omics dataset for a within -omics analysis (e.g. gene–gene, metabolite–metabolite).
  - two -omics datasets for paired -omics analysis (e.g., gene-miRNA, gene-metabolite).
- Parameters Used:
  - correlation method: Method to compute the correlation coefficient in `createVectors()`. Default method is "spearman". Other options are "pearson", "bwmc" and "sparcc".
  - components: Number of components for the fitted mixture model in `discordantRun()`. Default is 3. Another option is 5.
  - subsampling: Whether to subsample the correlation coefficients to be fitted in `discordantRun()` during the model fitting. Default is FALSE.

## Output:
No output file is generated. The output is a list of some elements:
- A matrix with posterior probabilities of differential co-expression for each gene pair under the main diagonal.
- A vector with the posterior probabilities of differential co-expression for each gene pair.
- A class matrix of the differential co-expression based on the highest posterior probability.
- A class vector of the differential co-expression based on the highest posterior probability under the main diagonal.
- probability matrix of posterior probabilities of different differential co-expression classes for every gene pair.
- The log likelihood of the model fitting.

***Note:*** A matrix and a vector could carry the same information. The matrix contains the posterior probabilities of differential co-expression for each gene pair under the main diagonal and NA above the main diagonal to eliminate redundancy. The vector is a vectorized (flattened) version of the matrix.

## Dependencies:
> - R version = **4.4**
> - BiocManager
## Commands used to install the tool:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("discordant")
```

## Commands used to run the tool:
- R code:
   - `createVectors()` function: Computes the correlation coefficients for every gene pair in both conditions.
   - `discordantRun()` function: Fits a mixture model, using the Expectation Maximization algorithm, to the correlation coefficients of every condition. Then, it estimates posterior probabilities as a significance test for differential co-expression.

## Reason the tool is impossible to install:
- None was found.

