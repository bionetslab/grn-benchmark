# pDNA

pDNA algorithm (https://doi.org/10.1093/bioinformatics/btx208), implemented in Matlab on github (https://github.com/Zhangxf-ccnu/pDNA), Using nonparanormal graphical models to analyze differential networks, the pDNA tool provides Matlab code. With this method, gene expression data can be analyzed more effectively, providing a robust framework for the identification of differential networks.


## Description

An approach to identify differential networks in biological research that is based on nonparanormal graphical models can be used by adding historical data to the interpretation of gene expression data. Adding historical data improves the interpretation of gene expression data and provides a solid foundation for finding differential networks.

This method overcomes the limitations of traditional differential network analysis by integrating prior knowledge, which allows for a more accurate analysis of network changes across different environments. This approach is particularly useful for analyzing high-dimensional, noisy biological data.

Sample Nonparanormal Covariance Matrices: These matrices are computed by the program to make non-Gaussian data analysis easier.
pDNA Model Solution: The tool solves the pDNA model.
Demos: Includes practical demonstrations using ovarian cancer and glioblastoma gene expression data to help users test the algorithm's effectiveness.

## Requirements

MatlabR2024a


## Installation Instructions

install MatalbR2024a

clone the pDNA tool

clone the data sets
```
git clone https://github.com/bionetslab/grn-benchmark
cd src/pDNA-tool 
```


## Instructions:

1. Open MATLAB.
2. Navigate to the directory containing pDNAtool.m.
3. Run the script by typing pDNAtool in the MATLAB command window.
4. Follow the prompts to input the directories for the two data files and the lambda value.
Example 
```
>> pDNAtool
Please put the first directory: path_to_first_file.tsv
Please put the second directory: path_to_second_file.tsv
Please choose your lambda: 0.45


```

## Relevant Parameters:
Lambda : You can experiment with different lambda values to see how it affects the sparsity and structure of the differential network.
for the refrence datsetets which were applied to the tool the value of lambda after many times was estimated to 0.45


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
The script will output a file named network.tsv containing the significant gene interactions.
- `target`: Targeted genes
- `regulator`: Regulatory genes
- `condition`: Condition that the edge belongs to
- `weight`: Weight of the edge

It will also display a heatmap of the weighted differential network sorted by node degrees.

## Troubleshooting
- Ensure the input file paths are correct and the files are in TSV format.
- Ensure Sigma_compute and pDNA functions are accessible in the MATLAB path.
- adjust the lambda value accordingly.




This README.md file was generated on 22.05.2024 by Ahmed Elshabassy (ahmed.elshabassy@fau.de)
