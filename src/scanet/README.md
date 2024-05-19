# SCANet

SCANet is a Python package designed to infer gene co-expression networks, extends them to gene regulatory networks and identifies potential drug candidates from scRNA-seq data.

## Description

SCANet offers an integrated approach to analyze single-cell gene expression data. Initially, it partitions the data into subclusters to extract representative cells.Subsequently, utilizing Spearman, Pearson, or biweight correlations, it infers Gene Co-expression Networks (GCNs). These GCNs are further analyzed to extract modules of interest, from which Gene Regulatory Networks (GRNs) are inferred. Ultimately, SCANet identifies potential drug candidates.

[SCANet Paper](https://academic.oup.com/bioinformatics/article/39/11/btad644/7325353)

[SCANet Github](https://github.com/oubounyt/SCANet)


## Installation Instructions

### Alternative 1) Docker:

- Build Docker image using `Dockerfile` in the `Docker` folder :
  
```bash
docker build -t scanet-app .
```
- Run the image using the following command:

```bash
docker run --rm -it scanet-app bash
```

### Alternative 2) Conda:

- Install environment: 
```bash
conda env create -f scanet.yml
```
- Activate environment:
```bash
conda activate scanet
```
- Find the scanet folder using the following command :
  
```bash
pip show scanet
```

-  Navigate to scanet folder and replace the file `grn.py` with the modified one, which is named `grn.py` and can be found in the `Docker` folder. This is important as I have changed this file slightly (Further details provided in the `Problem with grn inference function of SCANet` section).

- Please refer to the log file if you encounter any issues while running the tool.

## Examplery Execution Instructions:

- Docker :

Once inside the Docker container, execute the respective Python scripts using the following commands:

For `scanet_on_ref_data.py`:

```bash
python scanet_on_ref_data.py --c "/app/data/reference_datasets/full_input/out_CD8_exhausted.tsv" --d "/app/data/reference_datasets/full_input/out_Macrophages.tsv" --o "output/" --n 75 --t 0.8
```

For `first_study.py`:

```bash
python first_study.py --d "/app/data/GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad" --c "/app/data/Supplementary file 3.csv" --o "output/"
```

For `second_study.py`:

```bash
python second_study.py --d "/app/data/GSE147319_adata_villus_processed.h5ad" --o "output/"
```
These commands will run the specified Python scripts with the provided arguments, generating output files in the `output/` directory.

- Without docker:

```bash
python scanet_on_ref_data.py --c "..\..\reference_datasets\full_input\out_CD8_exhausted.tsv" --d "..\..\reference_datasets\full_input\out_Macrophages.tsv" --o "output" --n 75 --t 0.8
```

Note: Other datasets (e.g., 500 genes, 1000 genes, 2500 genes) can also be used by specifying their file paths.

--n argument is the number of representative cells/samples and --t is the correlation cutoff. Both are not required as they have default values 75 and 0.8 respectively.

## SCANet parameters

- `N` : Number of representative cells/samples per cluster.
- `correlation_cutoff` : Specifies the minimum threshold for co-expression strength that defines whether an edge should be included in the network. Edges with co-expression values below this threshold will be excluded from the network.


## Input file format specification:

- `input file 1`:  Path to a tab-separated file that contains the normalized gene expression for condition 1:
    - Each row represents a gene, while each column represents a cell. 
    - The first column, labeled **Gene**, lists the gene names.
    - The entries within the table denote normalized gene expression values.
- `input file 2`: Path to tab-separated file that contains the gene expression dataframe for condition 2:
    - Each row represents a gene, while each column represents a cell. 
    - The first column, labeled **Gene**, lists the gene names.
    - The entries within the table denote normalized gene expression values.
- `output path`: String specifying the directory where the output will be saved. It's essential that this directory exists before running the script.

## Output file format specification:

Output is a tab-separated file `network.tsv` stored at `output path`. It contains the following columns:
- `target`: Target genes
- `regulator`: Regulatory genes
- `condition`: Condition that the edge belongs to
- `weight`: Weight of the edge, which is the correlation strength.


`net.pkl` : A pickle file that include the inferred modules. It will be used later in the downstream analysis.

## Interpretation of the output

- Nodes represent genes.
- In GRN plots, triangles denote transcription factors, while stars represent potential drug candidates. In GCN plots, red nodes indicate hub genes. The blue nodes indicate target genes.
- Each edge indicates a notable interaction between a transcription factor (regulator) and a target gene, determined by Spearman, Pearson, or biweight correlation methods.

## Recommended hyperparameters by the authors

- The author published a pdf file `Supplementary file 1` to offer a guideline, how to determine the optimal number of subclusters, the impact of 'N' on Analysis Outcome and how to identify the modules of interest. I uploaded the file here in this repo.

- The defualt value of the `correlation_cutoff` is 0.6.

## Problem with grn inference function of SCANet

- The function grn_inference is responsible for generating the grn_df for a specified module. It requires several arguments:

sn.grn.grn_inference(adata_processed=adata_processed,modules_df=modules_df, module=Mod_, groupby_=groupby_, anno_name=anno_name, specie=specie_, subsampling_pct=80, n_iteration=n_iteration, num_workers=num_workers)

- Despite trying various parameter settings and on different pcs, the function consistently reports that 0 edges were found. Even when excuting the provided example with the same parameters, the result remains the same.

- After some debugging, I found out that the problem is in the grn.py file, line 262. There is a problem with the multiprocessing library so the function `regulon` is never excuted. I tried to cancel the multiprocessing calling and make it nonparallel but it didn't work as the multiprocessing library is used in other multiple places in the code.

- I have created an issue on the SCANet GitHub repository, accessible [here](https://github.com/oubounyt/SCANet/issues/8), and am waiting for assistance. I will contact you as soon as the matter is resolved, even if the response comes after submission.

- To deal with this problem, I created the grn_df my self. Since the function grn_inference identifies without any issues the transcription factors within a given module, I modified the script slightly to return these transcription factors. Then, utilizing the GCN of the same module, I could get the target genes for each transcription factor. This approach allowed me to construct the grn dataframe. The modified grn.py file is uplaoded on this repo with the name `grn.py` and can be found in the `Docker` folder.