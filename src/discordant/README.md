# Discordant

Discodant algorithm (https://doi.org/10.1186/s13104-016-2331-9), implemented in the Discordant R package (https://bioconductor.org/packages/release/bioc/html/discordant.html), is an algoritm designed for detecting pairs of genes that show differential co-expression between two conditions.


## Description
The discodant algorithm computes correlation coefficients for every gene pair in both conditions and use Bayesian statistics to estimate posterior probabilities as a signicance test for differential co-expression.

It fits a mixture model, using the Expectation Maximization algrithm, to the correlation coefficients of every condition. The algorithm doesn't produce a network, but rather a list of gene pairs and their associated posterior probability of being differentially co-expressed between the two conditions. 

The output of the algorithm is converted into a network by creating an edge, based on correlation coefficients difference, between every gene pair that show significant differential co-expression.

## Installation Instructions

### Alternative 1. R:
Install the required R packages:
```
install.packages(c("tidyverse","data.table","optparse"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("discordant", "clusterProfiler", "org.Hs.eg.db"))
```

### Alternative 2. Conda:
Create conda environment: 
```
conda env create -f environment.yml
```
Activate the environment:
```
conda activate discordant
```
### Alternative 3. Docker:
Build Docker image:
```
docker build -t discordant -f Dockerfile .
```
Run Docker container:
```
docker run --rm -it discordant bash
```

## Examplery Execution Instructions:
- Without docker:
```
# Full Input
Rscript discodant.R -a ../../reference_datasets/full_input/out_CD8_exhausted.tsv -b ../../reference_datasets/full_input/out_Macrophages.tsv -o ./
```
**Notes**:
- The are muliple optional arguments that can be used to specify the parameters in the script with default values. For more details, run the following command: 
```
 Rscript discodant.R --help 
``` 


- The algorithm can run on any dataset with the specified input format. The output will be stored in the specified output folder.

## Relevant Parameters:
- `cor.method`: Method to compute the correlation coefficient in `createVectors()` function. Default method is "spearman". other options are "pearson", "bwmc" and "sparcc".
- `components`: Number of components for the mixture model in `discordantRun()` function. Default is 3. Another option is 5.
- `subsampling`: whether to subsample the correlation coefficients to be fitted in `discordantRun()` function during the model fitting. Default is FALSE.
- `iter`: Number of iterations for subsampling in `discordantRun()` function. Default is 100.
- `transform`: Whether to Fisher z-transform the correlation coefficients for the mixture model to fit in `discordantRun()` function. Default is TRUE.
- `pp_threshold`: Posterior probability threshold for significant differential co-expression. Default is 0.95.
- `abs_weight_thresh`: Absolute weight,correlation difference across the conditions, threshold for significant differential co-expression. Default is 0.6.


## Input file format specification:
- `--input_1` or `-a`: Path to tab-separated file that contains the gene expression dataframe for condition 1:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value
- `--input_2` or `-b`: Path to tab-separated file that contains the gene expression dataframe for condition 2:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value
- `--output_path` or `-o`: String that contains the path to the output folder. Has to exist at time of execution.

## Output file format specification:
Tab-separated output file `network.tsv` stored at `--output_path` containing the following columns:
- `target`: Targeted genes
- `regulator`: Regulatory genes
- `condition`: Condition that the edge belongs to
- `weight`: Weight of the edge

## Interpretation of the output
- The nodes correspond to the genes.
- Each edge represents an interaction between a pair of genes that show differential co-expression(according to the discordance posterior probabilities) in either condition 1 or 2.

## Recommended hyperparameters by the authors
- `cor.method = "spearman"`, Spearman's correlation has the best performance.
- `components = 3`, 3 components are simpler yet sufficient to model the correlation coefficients and has better performance and faster runtime compared to 5 components.

## Neccessary information for the execution of the algorithm
Running this algorithm on a device with Amd Ryzen 5 5600h 6-Core Processor 3.59 GHz and 16 GB of RAM got the following runtimes for the different reference datasets:

| reference_dataset | subsampling | runtime(in minutes) | get_output |
| :---------------: | :---------: | :-----------------: | :--------: |
|     500 genes     |    False    |          4          |    True    |
|     500 genes     |    True     |          1          |    True    |
|    1000 genes     |    False    |         15          |    True    |
|    1000 genes     |    True     |          2          |    True    |
|    2500 genes     |    False    |         50          |   False    |
|    2500 genes     |    True     |          8          |    True    |
|    full input     |    False    |         80          |   False    |
|    full input     |    True     |         13          |    True    |

For the large datasets using the mentioned specs, the algorithm failed to generate output, due to the high computational cost, unless subsampling was used leading to reduce computational cost and runtime while preserving similar performance. Thus, `subsampling` is recommended for large datasets and is set as `True` in **discodant.R** script.

This README.md file was generated on 08.05.2024 by Karim Sakr (karim.m.sakr@fau.de)