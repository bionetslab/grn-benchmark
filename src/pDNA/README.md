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


## Examplery Execution Instructions:


## Relevant Parameters:


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



This README.md file was generated on 22.05.2024 by Ahmed Elshabassy (ahmed.elshabassy@fau.de)