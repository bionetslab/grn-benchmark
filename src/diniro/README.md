# DiNiro

DiNiro is a web based tool used in Differential co-expression regulatory networks analysis for Single-cell

## Description

DiNiro is developed for conjoint de novo differential scGRN reconstruction and network enrichment aims to regulate molecular mechanisms directly from scRNAseq data. DiNiro is able to get variation in co-expression patterns across two conditions based on gene expression profiles.

[DiNiro Tool] (https://apps.cosy.bio/diniro/)

[DiNiro Paper](https://academic.oup.com/nargab/article/5/1/lqad018/7069287)

[SCANet Github](https://gitlab.com/mhanedd/diniro)


## Installation Instructions

1) Docker:

- You can build Docker image using `Dockerfile` in the `Docker` folder :
  
```bash
docker compose build
```

- Then run the image using the following command:

```bash
docker compose up
```

- Access the tool through the local host link http://localhost:8025.


## DiNiro Usage

1. First you need to transfer the reference dataset to .H5AD file format using the script `Transfer Data to H5AD Script.ipynb` 

2. Open DiNro local host: http://localhost:8025 after you compose up the docker file.

3. Click on Get Started to go to the upload data.

4. Then Upload your file in `Anndata H5AD file`
5. If you have a file with specific transcription factors you can upload it in `TFs file` but it's optional
6. Add a name to the `Session Name` and click `Upload` to upload the data.

7. Choose `Plot Type` (Dimensionality Reduction Methods used to visualize the data like PCA and UMAP)
8. Choose `Plot Color` (Coloring Criteria to color cells)
9. Click on `Update Plot` to update it to the choosen type and color and then click `Next`. 

10. Select the sample that you would like to analyze either by `Cluster Selection` or `Lasso Selection` and click Next.

11. Add you parameters or use the default ones and click start to start the computations.

12. You will see a token link and a progress bar. After the progress bar reaches 100 %. Click on the access token link to see the results page.


## Examplery Execution Instructions(Background Tasks):

- After applying all the steps in the Usage section, let's see how the computations are done.

- Starts by retrieving matrices related to gene regulatory networks, which include interaction matrices between transcription factors and target genes based on the parameters provided 

- Then create `/data_folder_Sample/` for both samples to contain outputs.

```bash
        ComputeGrn.get_matrices(userrr, n_sub, s_size, TFs_type)
        current_user = userrr

        path_sample_1 = os.getcwd() + '/media/' + str(current_user) + '/data_folder_sample_1/'
        path_sample_2 = os.getcwd() + '/media/' + str(current_user) + '/data_folder_sample_2/'
```

- Then applies some data filtiring and processings and create another folder to contain the whole data that will be generated.

```bash
        sample_1 = FilteringProcessing.merging_and_filtring(path_sample_1, occurrence_threshold, size_thresholds)
        sample_2 = FilteringProcessing.merging_and_filtring(path_sample_2, occurrence_threshold, size_thresholds)
        data_folder_all = os.getcwd() + '/media/' + str(current_user)
```

- Then extracts these 2 files `sample_grn_s1.csv` `sample_grn_s2.csv` that will be used to continue the computations.

```bash
        sample_1.to_csv(path.join(data_folder_all, 'sample_grn_s1.csv'), index=False)
        sample_2.to_csv(path.join(data_folder_all, 'sample_grn_s2.csv'), index=False)
```
- Then feed the files extracted to Function.One to test data and extract regulons.

```bash

FunctionOne.functionOne(data_folder_all,size_thresholds)

```

- Then starts to apply the copula function to the data after analyzing relationships between transcription factors and other variables using empirical copulas

```bash

FunctionTwo.copula(data_folder_all)

```

- P-values and cutoff filtering

```bash
    p_values = Pvalues.p_values(data_folder_all)
```

- Then applying GeneEdgesDist to filters interaction data based on a p-value cutoff and handles file operations.

```bash
    GeneEdgesDist.update(data_folder_all,pval_cutoff)
```

- Lastly it starts to plot the networks using the all files required from `/data_folder_all/` that was created before
```bash

        Nlist, Elist, nlist, elist, n_combined, e_combined = Plotting.plot_network(data_folder_all, number_of_modules,pval_cutoff,int(size_threshold_min),int(size_threshold_max))
```

## DiNiro parameters

- `Significance cutoff` : We use a permutation test to calculate the p-value for every gene-gene interaction found in the results. Here you can specify the cut-off value to reject the null hypothesis, Default = 0.05.
- `Species selection` : We use transcription factors (TFs) when searching for gene regulation. In case you have not uploaded your own TFs of interest, species-specific transcription factor files are used (Here you can specify which file to use).
- `Number of subsamples` : Number of subsamples to derive from each sample in order to generate multiple gene regulatory networks (for denoising purpose), Default = 4
- `Sub-sampling size (%)` :The percentage of cells that will be random subsample from each sample to form the subsamples, Default = 70 %.
- `Number of Modules` : Number of gene modules to keep in the final results. Modules are ignored based on ranking criteria. 
- `Max module size` : The maximum number of genes in a module in order to include it in the results.
- `Occurrence threshold (%)` : To reduce noise, gene regulatory networks are computed multiple times and gene interactions are filtered based on an occurrence threshold (in percentage), Default = 70%.
- `Min module size` : The minimum number of genes in a module in order to include it in the results.


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

## Interpretation of the output

- Nodes represent genes.
- In GRNs plots, Blue triangles represent transcription factors.

## Recommended hyperparameters by the authors

- The author of the paper published a pdf file `Supplementary file 1` have some guidelines about choosing parameters. It will be provided in the repo in `Additional Files`

## Problems with the code
- Errors faced are mentioned in log.txt
- I have contacted the author of the tool for the latest source code which works on the live version of the website but I didn't have a reply yet.
