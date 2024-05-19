# SCANet Paper results Reproduction

This repository contains scripts for reproducing results of the two studies of the [SCANet](https://academic.oup.com/bioinformatics/article/39/11/btad644/7325353) paper.
The author hasn't puplished the parmeters used to generate the results of the paper. The inferred modules were published instead but this is not enough to reproduce the results. I contacted him to ask for the parmeters but got no reply. So it wasn't possible to reproduce exactly the same results. 
I worked on the same data that was used in the paper, followed the same methodology and workflow of SCANet, used the defult parameters and generated the results.

## Dependencies

- Python 3.6 or higher
- scanet
- pandas
- matplotlib
- os
- argparse

## Setup

- To run both scripts, ensure that the necessary Python packages are installed.
- Download the data used for both studies.

## Data

- [Cillo covid19 processed data GSE180578](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180578) was used in the first study. You can download it directly [here](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180578/suppl/GSE180578%5Fcillo%5Fcovid19%5Fstudy%5Faggregrated%5Fannotated%5Fdata%2Eh5ad%2Egz)

- [Villus processed data GSE147319](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147319) was used in the second study. You can download it directly [here](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147319/suppl/GSE147319%5Fadata%5Fvillus%5Fprocessed%2Eh5ad%2Egz)

## Scripts

### First Study (`first_study.py`)

This script processes `Cillo covid19 processed data`. It requires three input files: an h5ad data file, a CSV file containing modules data provided by the author, which can be find in this repo with the name `Supplementary file 3`, and an output directory where the output will be saved.

- **Usage**:
  ```bash
  python first_study.py --d <path_to_h5ad_file> --c <path_to_modules_csv> --o <output_directory>
  ```
  **Example**:
  ```bash
  python first_study.py --d "..\..\GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad" --c "..\Supplementary file 3.csv" --o ".\"
  ```

### Second Study (`second_study.py`)
This script processes `Villus processed data GSE147319`. It requires two input files: an h5ad data file and an output directory where the output will be saved.
- **Usage**:
  ```bash
  python second_study.py --d <path_to_h5ad_file> --o <output_directory>
  ```
  **Example**:
  ```bash
  python second_study.py --d "..\..\GSE147319_adata_villus_processed.h5ad" --o ".\"
  ```

## Output

In the provided output directory you will find the tsv files of the inferred grns and a folder called `outs_scan`, which contains the following plots :

**For first study**:
- `gcn_network_deceased` : The gcn plot of the module that exhibits the highest correlation with the "deceased" condition.
- `gcn_network_survived` : The gcn plot of the module that exhibits the highest correlation with the "survived" condition.
- `GRN_net_deceased` : The grn plot of the module that exhibits the highest correlation with the "deceased" condition.
- `GRN_net_survived` : The grn plot of the module that exhibits the highest correlation with the "survived" condition.

**For second study**:
`gcn+drugs_HFD_ISC` : The gcn plot of the module that exhibits the highest correlation with the "HFD_ISC" cell type plus potential drug candidates.
`GRN_net_HFD_ISC` : The grn plot of the module that exhibits the highest correlation with the "HFD_ISC" cell type.

## Output Figures

Here are the output figures from the first study:

### 1. UMAP Clusters of Cell Types
![UMAP Clusters of Cell Types](figures\umap_clasters_cell_types.png)

### 2. Visualization of classical monocytes
![UMAP Visualization of Patient Outcomes](figures\classical_monocytes.png)

### 3. Bar Graph of Genes per Module
![Bar Graph of Genes per Module](figures\genes_per_module.png)

### 4. Modules-Trait Correlation
![Modules-Trait Correlation](figures\correlation.png)

### 5.  GCN of the gene module associated with deceased patients
![Gene Interaction Network](figures\gcn_network_deceased.png)

### 6.  GCN of the gene module associated with survived patients
![Gene Interaction Network](figures\gcn_network_survived.png)

### 7.  GRN of the gene module associated with deceased patients
![Gene Interaction Network](figures\GRN_net_deceased.png)

### 8.  GRN of the gene module associated with survived patients
![Gene Interaction Network](figures\GRN_net_survived.png)



Here are the output figures from the second study:

### 1. UMAP Clusters
![UMAP Clusters](figures\clusters.png)

### 2. Modules-Trait Correlation
![Modules-Trait Correlation](figures\mod_trait_correlation.png)

### 3. GCN of the gene module associated with HFD_ISC + Drug candidates
![GCN of the gene module associated with HFD_ISC + Drug candidates](figures\GCN_net_HFD_ISC+drug.png)

### 4. GRN of the gene module associated with HFD_ISC
<<<<<<< HEAD
![GRN of the gene module associated with HFD_ISC](figures\GRN_net_HFD_ISC.png)
=======
![GRN of the gene module associated with HFD_ISC](GRN_net_HFD_ISC.png)
>>>>>>> 49e22c49f819fc2a9e744594643971ba4a63467f
