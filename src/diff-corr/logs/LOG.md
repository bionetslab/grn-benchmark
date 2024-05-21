# Tool Name: Differential Correlation Analysis using DiffCorr in R

## Function URL
The `DiffCorr` package can be found on CRAN or Bioconductor repositories.

## Input
- **input_file1**: File path to the first gene expression dataset (tab-separated values, TSV format).
- **input_file2**: File path to the second gene expression dataset (tab-separated values, TSV format).
- **output_file**: Desired string for the TXT file that will be generated, containing the differential correlation results.

## Output
- **DiffCorr_res.txt**: A text file containing the results of the differential correlation analysis.

## Dependencies
- `affy`
- `genefilter`
- `igraph`
- `spatstat`
- `cluster`
- `DiffCorr`

## Installation Commands
To install the required R packages, use the following commands in your R environment:

```R
install.packages("affy")
install.packages("genefilter")
install.packages("igraph")
install.packages("spatstat")
install.packages("cluster")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffCorr")
```

## Run Commands
```R
Rscript diffcorr_analysis.R input_file1.tsv input_file2.tsv DiffCorr_res.txt
```
