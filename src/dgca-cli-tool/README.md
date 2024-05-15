# Differential Gene Correlation Analysis Command-line Tool (DGCA CLI Tool)

<details>
<summary>Table of Contents</summary>

- [Brief Description](#brief-description)
- [Reference](#reference)
- [Available Commands](#available-commands)
- [Installation Instructions](#installation-instructions)
  - [Using Docker](#using-docker)
  - [Using Locally](#using-locally)
- [Parameters](#parameters)
- [Input File Format](#input-file-format)
- [Output File Format](#output-file-format)

</details>

## Brief Description
This tool is a command-line utility designed as a wrapper around the DGCA package, facilitating Differential Gene Correlation Analysis between two conditions. It aids in identifying significant changes in gene correlation across different biological conditions, supporting advanced biological insights into gene regulatory mechanisms.

## Reference
For detailed methodology and application, please refer to:
McKenzie, A.T., Katsyv, I., Song, WM. et al. DGCA: A comprehensive R package for Differential Gene Correlation Analysis. BMC Syst Biol 10, 106 (2016). [https://doi.org/10.1186/s12918-016-0349-1](https://doi.org/10.1186/s12918-016-0349-1).

## Available Commands
The CLI includes commands for:
- Differential Correlation Analysis (dgca.R)
- [Enrichment Map (enrichment_map.R)](downstream_analysis/enrichment-map.md)

In this readme, we will see the usage of dgca.R

## Installation Instructions

### Using Docker
To run the tool using Docker, ensure Docker is installed on your system and follow these steps:

1. Navigate to the project (dgca-tool) root directory.

2. Build the Docker image:
   ```bash
   docker build -t dgca-tool .
   ```
3. Run the tool using the Docker container:
   ```bash
   docker run --rm -v ./data:/data dgca-tool dgca.R --input_file_1 /data/BRCA_normal.tsv --input_file_2 /data/BRCA_tumor.tsv --output_path /data
   ```

### Using Locally (Using R)
1. Navigate to the project (dgca-tool) root directory.

2. Install the required R packages:
   ```bash
   Rscript install_dependencies.R
   ```
3. Run the tool:
   ```bash
   Rscript dgca.R --input_file_1 /data/BRCA_normal.tsv --input_file_2 /data/BRCA_tumor.tsv --output_path ./data
   ```

## Parameters
- `--input_file_1`: Path to the TSV file containing gene expression data for the first condition.
- `--input_file_2`: Path to the TSV file containing gene expression data for the second condition.
- `--output_path`: Output path where the results will be saved. The results file will be named 'network.tsv'.

## Input File Format
The input files should be TSV (Tab-Separated Values) format containing gene expression data. Each file must have genes as rows and samples as columns for each condition.

## Output File Format

| Target  | Regulator | Condition    | Weight_1      | Weight_2   | Weight_3     | Weight_4     | Weight_5 |
|---------|-----------|--------------|---------------|------------|--------------|--------------|----------|
| CACYBP  | NACA      | condition1   | -0.070261455  | 0.67509118 | 1.100991e-24 | 1.100991e-24 | 0/+      |
| CACYBP  | NACA      | condition2   | 0.9567267     | 0          | 1.100991e-24 | 1.100991e-24 | 0/+      |

### Column Descriptions:

- **Target**: Target node of the edge.
- **Regulator**: Source node of the edge.
- **Weight\_1 (Correlation Coefficients)**: Quantifies the strength and direction of the linear relationship between two genes within specified conditions. Depending on the distribution assumption, either the Pearson product-moment correlation coefficient ($r_p$) or the Spearman’s rank correlation coefficient ($r_s$) is used.
  
- **Weight\_2 (P-values for Correlations)**: Determines the statistical significance of the correlation coefficients, indicating the likelihood of observing the calculated correlations by chance under the null hypothesis of no association.
  
- **Weight\_3 (Z-Score Differences)**: Emphasizes the extent of variation in the correlation between two gene expressions across different conditions. It quantifies how much the relationship between these genes changes from one condition to another. The z-score difference is derived using the Fisher z-transformation, which normalizes the variance of correlation coefficients:
  
  $z = \frac{1}{2} \log_e \left(\frac{1+r}{1-r}\right)$

  where $r$ is the correlation coefficient, and $\log_e$  is the natural logarithm. The z-score difference is calculated as:
  
  $dz = \frac{(z_1 - z_2)}{\sqrt{|s_{z_1}^2 - s_{z_2}^2|}}$
  
  where $s_{z_x}^2$ represents the variance of the z-score in condition $x$.
  
- **Weight\_4 (P-values of the Z-Score Differences)**: Assesses the significance of the differences between the z-scores of two conditions using p-values.
  
- **Weight\_5 (Classifications of the Correlation Change)**: Evaluates whether gene pairs experience a gain or loss of correlation between two conditions. This classification further categorizes these changes based on their statistical significance and the directionality (positive or negative) of the correlation.
