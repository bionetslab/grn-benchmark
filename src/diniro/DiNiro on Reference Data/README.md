# DiNiro on the reference datasets

To feed the reference dataset you need to run DiNiro in the localhost and then it will be easy to deal with the local web interface.


## Requirements

- Reference dataset should be Anndata (.H5AD) file format.
- Transcription Factor file in (.txt) file format


## Method parameters

- `Significance cutoff` : We use a permutation test to calculate the p-value for every gene-gene interaction found in the results. Here you can specify the cut-off value to reject the null hypothesis, Default = 0.05.
- `Species selection` : We use transcription factors (TFs) when searching for gene regulation. In case you have not uploaded your own TFs of interest, species-specific transcription factor files are used (Here you can specify which file to use).
- `Number of subsamples` : Number of subsamples to derive from each sample in order to generate multiple gene regulatory networks (for denoising purpose), Default = 4
- `Sub-sampling size (%)` :The percentage of cells that will be random subsample from each sample to form the subsamples, Default = 70 %.
- `Number of Modules` : Number of gene modules to keep in the final results. Modules are ignored based on ranking criteria. 
- `Max module size` : The maximum number of genes in a module in order to include it in the results.
- `Occurrence threshold (%)` : To reduce noise, gene regulatory networks are computed multiple times and gene interactions are filtered based on an occurrence threshold (in percentage), Default = 70%.
- `Min module size` : The minimum number of genes in a module in order to include it in the results.

## Usage

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

13. Results Page have different options to choose from:
    - `Table`: the table have the following columns: Rank, Module TF, Genes, GO Terms, Module Network
    - `Network`: Visualize all modules together and you can see the genes regulated by more than TF.
    - `Summary`: Where you can download a summary for your results as .CSV file that have distribution of the differentially expressed genes across the samples.


## Inputs

- `Anndata H5AD file`: your reference data in .H5AD file format.
- `TFs file`: Transcription Factor file in (.txt) file format.

## Outputs

- In thedirectory: \DiNiro\Docker\media\CrEJEGgRfYYaoB0FIsfgkPCi4Ew2Yet5
- CrEJEGgRfYYaoB0FIsfgkPCi4Ew2Yet5 --> random user created each time you run the tool has all the generated data from the tool.
- You will find the following : sample_grn_s1.csv, sample_grn_s2.csv

- Use the script `Script to Produce the Output.ipynb` to transform it to the required output shape then you will get the following files:

    - `network.tsv` : A tsv file which represents the two GRNs of the two module, which have the highest correlation to macrophages and exhausted conditions.
    - `additional_data.tsv` : A tsv file have more additional data generated from the tool.
