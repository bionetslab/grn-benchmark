This README.md file was generated on 26/3/2024 by Nicolai Meyerhoefer (nicolai.meyerhoefer@fau.de)
# Zscore

Z-score method (https://doi.org/10.1093/bioinformatics/btm482) implemented via dcanr (https://www.bioconductor.org/packages/release/bioc/html/dcanr.html) is a method for inference of differential gene regulatory networks between two conditions.

## Description
The z-score differential gene regulatory network inference tool computes Pearson's or Spearman's correlations for each pair of genes for both conditions and applies Fisher z-transformation afterwards. The resulting z-scores are modeled as a normal distribution and a z-test is used to detect significant pairwise edges.

## Installation Instructions
### Alternative 1) Conda:
Install environment: 
```
conda env create -f conda.yml
```
Activate environment:
```
conda activate conda.yml
```
### Alternative 2) Docker:
Build Docker image:
```
docker build -t zscore -f Dockerfile .
```

## Examplery Execution Instructions:
- Without docker:
```
# 500 genes
Rscript zscore.R -c ../../data/reference_datasets/500_genes/out_CD8_exhausted.tsv -d ../../data/reference_datasets/500_genes/out_Macrophages.tsv -o ./

# 1000 genes
Rscript zscore.R -c ../../reference_datasets/1000_genes/out_CD8_exhausted.tsv -d ../../reference_datasets/1000_genes/out_Macrophages.tsv -o ./

# 2500 genes
Rscript zscore.R -c ../../reference_datasets/2500_genes/out_CD8_exhausted.tsv -d ../../reference_datasets/2500_genes/out_Macrophages.tsv -o ./

# Full Input
Rscript zscore.R -c ../../reference_datasets/full_input/out_CD8_exhausted.tsv -d ../../reference_datasets/full_input/out_Macrophages.tsv -o ./
```

## Input file format specification:
- `--input.file.1`: Path to tab-separated file that contains the gene expression dataframe for condition 1:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value
- `--input.file.2`: Path to tab-separated file that contains the gene expression dataframe for condition 2:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value
- `--output.path`: String that contains the path to the output folder. Has to exist at time of execution.

## Output file format specification:
Tab-separated output file stored at `--output.path` containing the following columns:
- `target`: Targeted genes
- `regulator`: Regulatory genes
- `condition`: Condition that the edge belongs to
- `weight`: Weight of the interaction

## Interpretation of the output
- The nodes correspond to the genes.
- Each edge represents an interaction between a pair of genes that is significant (according to the z-test) in either condition 1 or 2.
