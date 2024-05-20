# DCloc and DCglob

Identification of differentially correlated genes via local or global topology analysis (https://doi.org/10.1186/1752-0509-7-78). 


## Description

The dcloc_dcglob function applies a local or global topology analysis on the two given conditions and outputs a list of differentially correlated genes. These are then processed and combined to compute the differential gene regulatory network.

### 1) DCglob
Dcglob performs a global topology analysis by comparing connected components of the correlation networks of the given conditions. Genes that are present in connected components for only one of the conditions' correlation networks are considered as differentially correlated. The significance of the differential correlation is described by a $p$-value. A lower $p$-value corresponds to a stronger differential correlation. To obtain a list of differentially correlated genes a threshold `--pthresh` has to be stated. The default value is set to $0.1$. 

### 2) DCloc

Dcloc performs a local topology analysis by comparing the neighborhoods of the individual genes. The significance of the differential correlation is described by the metric topological dissimilarty $d$, in which the common next neighbors of a specific gene in both networks are compared. A higher $d$-value corresponds to a stronger differential correlation. To obtain a list of differentially correlated genes a threshold `--dthresh` has to be stated. The default value is set to $0.25$. 

## Installation Instructions
### Alternative 1) Conda:
Install environment: 
```
conda env create -f environment.yml
```
Activate environment:
```
conda activate dcloc_dcglob
```
### Alternative 2) Docker:
Build Docker image:
```
docker build -t dcloc_dcglob .

```
Run Docker container in interactive mode
```
docker run -it dcloc_dcglob /bin/bash

```
Go to work directory 

```
cd /home/analysis
```

Before executing shell script add as executible

```
chmod +x example_script.sh
```


## Important Files

- `12918_2013_1140_MOESM1_ESM.r`: Original code of Bockmayr et al. (https://doi.org/10.1186/1752-0509-7-78), contains functions dcloc and dcglob

- `dcloc_dcglob.R`: script for generating the differential GRN based on dcloc/dcglob

- `dcloc_dcglob_runner.sh`: shell script example for running dcloc_dcglob.R 

- `downstream_analysis.R`: script for downstream analysis

- `downstream_analysis_runner.sh`: shell script example for running downstream_analysis.R

- `false_discovery_rate_estimation.R`: script for downstream analysis

- `false_discovery_rate_runner.sh`: shell script example for running downstream_analysis.R


## Local and Global topology Analysis / Script: dcloc_dcglob.R

### Exemplary Execution Instruction 
#### With conda environment:

- local topology analysis
```
Rscript dcloc_dcglob.R -c data/out_CD8_exhausted.tsv -d data/out_Macrophages.tsv -o output --local --dthresh 0.25 --corrthresh 0.5
```

- global topology analysis
```
Rscript dcloc_dcglob.R -c -c data/out_CD8_exhausted.tsv -d data/out_Macrophages.tsv -o output --global --pthresh 0.1 --corrthresh 0.5
```


###  Input file format specification:
- `--input.file.1`: Path to tab-separated file that contains the gene expression data frame for condition 1:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value
- `--input.file.2`: Path to tab-separated file that contains the gene expression data frame for condition 2:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value

- `--output.path`: String that contains the path to the output folder. Has to exist at time of execution.

- `--local`: flag to determine that the local topology analysis shall be conducted, per default FALSE

- `--dthresh`: threshold for the topological dissimilarity for the local topology dissimilarity, will only be used if local topology analysis also selected, default = 0.25

- `--global`: flag to determine that the global topology analysis shall be conducted, per default TRUE

- `--pthresh`: threshold for the p-value for the global topology dissimilarity, will only be used if global topology analysis also selected, default = 0.1

- `--corrthresh`: threshold for the accepted correlation between genes and when to add the edge to the correlation network, default = 0.5

### Output format specification:

Tab-separated output file 'network.tsv' stored at `--output.path` containing the following columns:
- `target`: Targeted genes
- `regulator`: Regulatory genes
- `condition`: Condition that the edge belongs to
- `weight`: Correlation of the genes
- `edge_type`: Signifies type of edge e.g. undirected

- for global topology analysis:
Tab-separated output file 'DCglob.tsv' stored at `--output.path` containing the output of the algorithm DCglob
contains the colomns: 
- `p-value`: p-value signifiying differetial correlation
- `High. Corr.`: either 1 or 2, signifiying the condition for which the gene shows a stronger correlation

- for local topology analysis:

Tab-separated output file 'DCglob.tsv' stored at `--output.path` containing the output of the algorithm DCloc

contains the colomns: 
- `abs.top.dissim`: absolute topological dissimilarity
- `top.dissim`: topological dissimilarity
- `mean.neigh.A`: mean neighborhood of condition 1
- `mean.neigh.B`: mean neighborhood of condition 2

### Interpretation of the output in the file netowrk.tsv
- The nodes correspond to the detected differentially correlated genes.
- Each edge represents the correlation between a pair of genes that are differentially correlated in either condition 1 or 2.


## Downstream Analysis / Script: downstream_anaylsis.R

This script conducts the downstream analysis based on the list of genes provided by the local and global topology analysis. The list of genes are sorted into two lists containing the genes with a stronger correlation for condition 1 and with genes with a stronger correlation for condition 2. These lists are then used to create networks, perform hierarchical clustering and plot heatmaps of the expression data. 

The processing of the data and creation of the plots follows the procedure introduced by Bockmayr et al..

The output contains:
- networks in '.gml'- format of the differentially correlated genes for both conditions 
- Hierarchical Clustering for both conditions
- Heatmaps with color bars visualizing the hierarchical clustering for both conditions


### Exemplary Execution Instruction 
- With conda environment:

```
Rscript downstream_analysis.R-c data/out_CD8_exhausted.tsv -d data/out_Macrophages.tsv -o output --dthresh 0.21 --pthresh 0.1 --corrthresh 0.5 --dcloc.output dcloc_output.tsv --dcglob.output dcglob_output.tsv 
```

###  Input file format specification:
- `--input.file.1`: Path to tab-separated file that contains the gene expression data frame for condition 1:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value

- `--input.file.2`: Path to tab-separated file that contains the gene expression data frame for condition 2:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value

- `--output.path`: String that contains the path to the output folder. Has to exist at time of execution. 

- `--dcloc.output`: Path to tab-separated file that contains the output of the DCloc algorithm

- `--dcglob.output`: Path to tab-separated file that contains the output of the DCglob algorithm

- `--dthresh`: threshold for the topological dissimilarity for the local topology dissimilarity, default = 0.25

- `--pthresh`: threshold for the p-value for the global topology dissimilarity, default = 0.1

- `--corrthresh`: threshold for the accepted correlation between genes and when to add the edge to the correlation network, default = 0.5


### Output format specification:
All files are stored at `--output.path` 

Networks:
- `graph_DCloc_cond.name_cond1.name.gml`: correlation network of cond1.name containing the genes from the local topology analysis with a stronger correlation for condition cond.name

- `graph_DCglob_cond.name_cond1.name.gml`: correlation network of cond1.name containing the genes from the global topology analysis with a stronger correlation for condition cond.name

Heatmaps:
- `Heatmap_DCloc_cond.name_cond1.name.pdf`: Heatmap of expression data of cond1.name containing the genes from the local topology analysis with a stronger correlation for condition cond.name

- `Heatmap_DCglob_cond.name_cond1.name.pdf`: Heatmap of expression data of cond1.name containing the genes from the global topology analysis with a stronger correlation for condition cond.name

Hierarchical clustering:
- `clusters_global_cond.name.txt`: contains the hierarchical clustering  for the genes from the global topology analysis with a stronger correlation for condition cond.name 
- `clusters_local_cond.name.txt`: contains the hierarchical clustering  for the genes from the local topology analysis with a stronger correlation for condition cond.name 


## False Discovery Rate Estimation / Script: false_discovery_rate_estimation.R

For the selection of the threshold dthresh and pthresh to detect
the differentially correlated genes of the local and
global topology analysis, a statistical evaluation by
a repeated random subsampling analysis was introduced by Bockmayr et. al. The script conducts the false discovery rate estimation for the provided thresholds and outputs the values in the terminal. 

### Exemplary Execution Instruction 

#### With conda environment:
```
Rscript false_discovery_rate_estimation.R -c data/out_CD8_exhausted.tsv -d data/out_Macrophages.tsv -o output  --dthresh 0.25 --pthresh 0.1 --n.samples 100 --n.repetitions 100
```

###  Input file format specification:
- `--input.file.1`: Path to tab-separated file that contains the gene expression data frame for condition 1:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value
- `--input.file.2`: Path to tab-separated file that contains the gene expression data frame for condition 2:
    - Rows correspond to gene names and columns to cells 
    - First column is named **Gene** and contains the gene names
    - Each entry is a normalized gene expression value

- `--output.path`: String that contains the path to the output folder. Has to exist at time of execution.
- `--dthresh`: threshold for the topological dissimilarity for the local topology dissimilarity, will only be used if local topology analysis also selected, default = 0.25
- `--pthresh`: threshold for the p-value for the global topology dissimilarity, will only be used if global topology analysis also selected, default = 0.1
- `--n.samples`: number of samples for the subsampled data, default = 100
- `--n.repetitions`: number of repetitions for n0 estimation, default = 100
