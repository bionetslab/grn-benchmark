# ChNet-Tool

# Differential network analysis (ChNet) Tool

## Brief Description
This tool is a command-line utility designed as a wrapper around the ChNet package that can do differential network analysis by simultaneously analyzing differential gene interaction and differential gene expression. This helps to detect notable variations in gene correlation between various biological states, providing a deeper understanding of gene regulatory processes.

## Installation Software
Download and install the software below
 - Install RStudio
 - Install Docker
   
## Installation Instructions
Clone and navigate to the project (ChNet) root directory.
   ```bash
   git clone git@github.com:bionetslab/grn-benchmark.git && cd grn-benchmark/src/ChNet-tool
   ```
You can run the project in two ways.
- Docker
- Rscript Command Line
  
### Using Docker
To run the tool using Docker, ensure Docker is installed on your system and follow these steps:

1. Navigate to the project (ChNet-tool) root directory.

2. Build the Docker image:
   ```bash
   docker build -t docker .
   ```
## Reference Data Input Output
### Execute the ChNet tool using Docker for different datasets

**500 Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" -v "/${PWD}/Ref_dataset:/usr/src/app/Ref_dataset" docker Rscript Reference_data_Input_output.R  --input1 "Ref_dataset/500/out_CD8_exhausted.tsv" --input2 "Ref_dataset/500/out_Macrophages.tsv" --output_path "output/TSV_files/500" --file_name_suffix "500"
```
**1000 Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" -v "/${PWD}/Ref_dataset:/usr/src/app/Ref_dataset" docker Rscript Reference_data_Input_output.R  --input1 "Ref_dataset/1000/out_CD8_exhausted.tsv" --input2 "Ref_dataset/1000/out_Macrophages.tsv" --output_path "output/TSV_files/1000" --file_name_suffix "1000"
```
**2500 Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" -v "/${PWD}/Ref_dataset:/usr/src/app/Ref_dataset" docker Rscript Reference_data_Input_output.R  --input1 "Ref_dataset/2500/out_CD8_exhausted.tsv" --input2 "Ref_dataset/2500/out_Macrophages.tsv" --output_path "output/TSV_files/2500" --file_name_suffix "2500"
```
**TCGA.BRCA Genes:**
```bash
 docker run --rm -v "/${PWD}/output:/usr/src/app/output" docker Rscript TCGA.BRCA_TSV.R --dataset "TCGA.BRCA" --output_path "output/TSV_files/paper"

```

**GSE13159.AML Genes:**
```bash
 docker run --rm -v "/${PWD}/output:/usr/src/app/output" docker Rscript GSE13159.AML_TSV.R --dataset "GSE13159.AML" --output_path "output/TSV_files/paper"
```

## Downstream Analysis for the ChNet tool
### Execute the ChNet tool using Docker for different datasets for downstream analysis

**500 Genes:**
```bash
docker run -v "/${PWD}/output:/usr/src/app/output" -v "/${PWD}/Ref_dataset:/usr/src/app/Ref_dataset" docker Rscript downstream_analysis.R --operation_type "ref" --input1 "Ref_dataset/500/out_CD8_exhausted.tsv" --input2 "Ref_dataset/500/out_Macrophages.tsv" --condition1 "CD8 Exhausted T-cells" --condition2 "Macrophages" --output_path "output/downstream/500" --pdf_name_suffix "500"
```
**1000 Genes:**
```bash
docker run -v "/${PWD}/output:/usr/src/app/output" -v "/${PWD}/Ref_dataset:/usr/src/app/Ref_dataset" docker Rscript downstream_analysis.R --operation_type "ref" --input1 "Ref_dataset/1000/out_CD8_exhausted.tsv" --input2 "Ref_dataset/1000/out_Macrophages.tsv" --condition1 "CD8 Exhausted T-cells" --condition2 "Macrophages" --output_path "output/downstream/1000" --pdf_name_suffix "1000"
```
**2500 Genes:**
```bash
docker run -v "/${PWD}/output:/usr/src/app/output" -v "/${PWD}/Ref_dataset:/usr/src/app/Ref_dataset" docker Rscript downstream_analysis.R --operation_type "ref" --input1 "Ref_dataset/2500/out_CD8_exhausted.tsv" --input2 "Ref_dataset/2500/out_Macrophages.tsv" --condition1 "CD8 Exhausted T-cells" --condition2 "Macrophages" --output_path "output/downstream/2500" --pdf_name_suffix "2500"
```
**TCGA.BRCA Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" docker Rscript downstream_analysis.R --operation_type "paper" --input1 "TCGA.BRCA" --condition1 "Basal" --condition2 "LumA" --output_path "output/downstream/paper" --pdf_name_suffix "TCGA.BRCA"
```

**GSE13159.AML Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" docker Rscript downstream_analysis.R --operation_type "paper" --input1 "GSE13159.AML" --condition1 "cancer" --condition2 "normal" --output_path "output/downstream/paper" --pdf_name_suffix "GSE13159.AML"
```
## Replicate the Paper
### Execute the ChNet tool using Docker to replicate the figures of ChNet paper

**TCGA.BRCA Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" docker Rscript replicate_paper.R --input1 "TCGA.BRCA" --condition1 "Basal" --condition2 "LumA" --output_path "output/replicate_paper" --pdf_name_suffix "TCGA.BRCA"
```

**GSE13159.AML Genes:**
```bash
docker run --rm -v "/${PWD}/output:/usr/src/app/output" docker Rscript replicate_paper.R --input1 "GSE13159.AML" --condition1 "cancer" --condition2 "normal" --output_path "output/replicate_paper" --pdf_name_suffix "GSE13159.AML"
```


### Running Locally (Using R)
 According to the [chNet](https://github.com/Zhangxf-ccnu/chNet) R package on GitHub, it should be installed as below.
```R
# Step 1. Install the devtools package. Invoke R and then type
install.packages("devtools") 

# Step 2. Load the devtools package.
library("devtools") 

# Step 3. Install the chNet package from GitHub.
install_github("Zhangxf-ccnu/chNet", subdir="pkg")
```

## Execute the Chnet analysis using Rscript directly for reference gene and paper gene datasets as follows:
**Execute while staying in the R folder** 

**500 Genes:**
```bash
Rscript Reference_data_Input_output.R --input1 "ref_dataset/500/out_CD8_exhausted.tsv" --input2 "ref_dataset/500/out_Macrophages.tsv" --output_path "Results/TSV_files/500" --file_name_suffix "500"
```
**1000 Genes:**
```bash
Rscript Reference_data_Input_output.R --input1 "ref_dataset/1000/out_CD8_exhausted.tsv" --input2 "ref_dataset/1000/out_Macrophages.tsv" --output_path "Results/TSV_files/1000" --file_name_suffix "1000"
```
**2500 Genes:**
```bash
Rscript Reference_data_Input_output.R --input1 "ref_dataset/2500/out_CD8_exhausted.tsv" --input2 "ref_dataset/2500/out_Macrophages.tsv" --output_path "Results/TSV_files/2500" --file_name_suffix "2500"
```
**TCGA.BRCA Genes:**
```bash
 Rscript TCGA.BRCA_TSV.R --dataset "TCGA.BRCA" --output_path "Results/TSV_files/paper"
```

**GSE13159.AML Genes:**
```bash
 Rscript GSE13159.AML_TSV.R --dataset "GSE13159.AML" --output_path "Results/TSV_files/paper"
```


## Execute the ChNet analysis using Rscript directly for downstream analysis as follows:

**500 Genes:**
```bash
Rscript downstream_analysis.R --operation_type "ref" --input1 "ref_dataset/500/out_CD8_exhausted.tsv" --input2 "ref_dataset/500//out_Macrophages.tsv" --condition1 "CD8 Exhausted T-cells" --condition2 "Macrophages" --output_path "Results/downstream/500" --pdf_name_suffix "500"
```

**1000 Genes:**
```bash
Rscript downstream_analysis.R --operation_type "ref" --input1 "ref_dataset/1000/out_CD8_exhausted.tsv" --input2 "ref_dataset/1000//out_Macrophages.tsv" --condition1 "CD8 Exhausted T-cells" --condition2 "Macrophages" --output_path "Results/downstream/1000" --pdf_name_suffix "1000"
```
**2500 Genes:**
```bash
Rscript downstream_analysis.R --operation_type "ref" --input1 "ref_dataset/2500/out_CD8_exhausted.tsv" --input2 "ref_dataset/2500//out_Macrophages.tsv" --condition1 "CD8 Exhausted T-cells" --condition2 "Macrophages" --output_path "Results/downstream/2500" --pdf_name_suffix "2500"
```
**TCGA.BRCA Genes:**
```bash
Rscript downstream_analysis.R --operation_type "paper" --input1 "TCGA.BRCA" --condition1 "Basal" --condition2 "LumA" --output_path "Results/downstream/paper" --pdf_name_suffix "TCGA.BRCA"
```

**GSE13159.AML Genes:**
```bash
 Rscript downstream_analysis.R --operation_type "paper" --input1 "GSE13159.AML" --condition1 "cancer" --condition2 "normal" --output_path "Results/downstream/paper" --pdf_name_suffix "GSE13159.AML"
```
## Replicate the ChNet paper analysis using Rscript directly for as follows:

**TCGA.BRCA Genes:**
```bash
 Rscript replicate_paper.R --input1 "TCGA.BRCA" --condition1 "Basal" --condition2 "LumA" --output_path "Results/replicate_paper" --pdf_name_suffix "TCGA.BRCA"
```

**GSE13159.AML Genes:**
```bash
 Rscript replicate_paper.R --input1 "GSE13159.AML" --condition1 "cancer" --condition2 "normal" --output_path "Results/replicate_paper" --pdf_name_suffix "GSE13159.AML"
```


## Parameters
- `--input_file_1`: Path to the TSV file containing gene expression data for the first condition.
- `--input_file_2`: Path to the TSV file containing gene expression data for the second condition.
- `--output_path`: Output path where the results will be saved. The results file will be named 'network.tsv'.

## Input File Format
The input files should be TSV (Tab-Separated Values) format containing gene expression data. Each file must have genes as rows and samples as columns for each condition.

## Output File Format
| Target  | Regulator | Condition        | Weight     |
|---------|-----------|--------------    |------------|
| CLCA2  | BAMBI      | non-differential | 1          | 
| DAPL1  | SHC4       | differential     | 1          | 

**Column Descriptions:**
* `network.tsv`: tab-separated file that contains all edges (row-wise) with the following columns:
- **Target**: Target node of the edge.
- **Regulator**: Source node of the edge.
- **Weight**: Correlation coefficient of the gene pair for the specified condition.
- **Condition**: is it differential or non-differential

## Methodology

Algorithm: Workflow of the ChNet Tool

![Screenshot 2024-05-22 022941](https://github.com/RifatPiyal/ChNet-Tool/assets/51060047/8e54c05b-9755-4c61-ad7b-55057af372ed)

## Computational Steps for Differential Network Analysis
**Requirements:**
- Input File 1: Gene expression dataset for condition 1.
- Input File 2: Gene expression dataset for condition 2.
- Output Path: Directory to save the results.

**Procedure:**

1. **Data Preparation:**
- Load Datasets: Read the gene expression datasets from the specified input files into data frames.
- Verify Data Consistency: Ensure both datasets have the same genes in the same order by intersecting the gene lists.
- Combine Datasets: Combine the datasets from both conditions into a single matrix.
- Run Analysis: Perform the differential network analysis using the chNet package.



**condition:**
Indicates whether the interaction between the gene pairs is differential or non-differential. The condition is determined based on the differential expression status of the genes involved:

-differential: Both genes in the pair are differentially expressed.

-non-differential: At least one of the genes in the pair is not differentially expressed.

**Output** 
a network.tsv file detailing the differential coexpression network, specifying target, regulator, condition, and weight for each gene pair.
Genes and edges that exceed the defined thresholds ( 𝜆 𝑖 λ i ​ and 𝜆 𝑖 𝑗 λ ij ​ , respectively) are classified as differential (Diff); otherwise, they are classified as non-differential (Non-diff).

## Additional Material 
- [chNet's GitHub](https://github.com/Zhangxf-ccnu/chNet)
- [chNet's author doc](https://github.com/Zhangxf-ccnu/chNet/blob/master/chNet_2.0.0.pdf)
- [chNet's Paper](https://academic.oup.com/bioinformatics/article/37/23/4414/6318772)


<br>
## Reference

<a id="1">[1]</a>
Wickham, Hadley. (2016).
ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York.
ISBN: 978-3-319-24277-4 Add to Citavi project by ISBN.
ggplot2

<a id="2">[2]</a>
Carlson, Marc. (2023).
org.Hs.eg.db: Genome wide annotation for Human.
R package version 3.18.0.

<a id="3">[3]</a>
Tu, Jia-Juan. (2024).
chNet: Differential network analysis by simultaneously considering changes in gene interactions and gene expression.
R package version 2.0.0.

<a id="4">[4]</a>
Auguie, Baptiste. (2017).
gridExtra: Miscellaneous Functions for "Grid" Graphics.
R package version 2.3.
gridExtra

<a id="5">[5]</a>
Huber, Wolfgang; Carey, V.J.; Gentleman, R.; et al. (2015).
Orchestrating high-throughput genomic analysis with Bioconductor.
Nature Methods, 12, pages 115-121.
BiocGenerics

<a id="6">[6]</a>
Wickham, Hadley. (2007).
Reshaping Data with the reshape Package.
Journal of Statistical Software, 21(12), pages 1-20.
reshape2


<a id="7">[7]</a>
Csárdi, Gábor; Nepusz, Tamás. (2024).
igraph: Network Analysis and Visualization in R.
DOI: 10.5281/zenodo.7682609, Add to Citavi project by DOI R package version 2.0.3.
igraph

<a id="8">[8]</a>
Gu, Zuguang; Gu, Lei; Eils, Roland; Schlesner, Matthias; Brors, Benedikt. (2014).
circlize implements and enhances circular visualization in R.
Bioinformatics, 30(19), pages 2811-2812.

<a id="9">[9]</a>
Kolde, Raivo. (2019).
pheatmap: Pretty Heatmaps.
R package version 1.0.12.
pheatmap

