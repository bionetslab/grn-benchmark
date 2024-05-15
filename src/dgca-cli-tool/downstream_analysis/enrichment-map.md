# Enrichment Analysis 

## Description
This script performs Gene Set Enrichment Analysis (GSEA) on provided network data to identify significant pathways related to differential gene correlations. The script can generate enrichment maps for combined data or for specific conditions, facilitating deeper insights into the biological significance of the observed correlations.


## Usage (using docker)
If you have already built the docker image for the dgca-tool outlined in the main readme.md, just run the following command:
```bash
    docker run --rm -v <path_to_your_data>:/data dgca-tool enrichment_map.R --network_file /data/network.tsv --output_path /data
```
Replace <path_to_your_data> with the directory path where network.tsv exists.

## Usage (without docker)
Run the script from the command line by providing the necessary options. Make sure that all the requirements in `install_dependencies.R` file are installed beforehand.
 Below is the basic usage format:
```bash
Rscript enrichment_map.R --network_file [path to TSV file] --output_path [path to save output]
```

### Options
- `-f`, `--network_file`: Path to the TSV file containing network data.
- `-o`, `--output_path`: Output directory to save the enrichment map image.
- `--ont`: Ontology to use for GSEA (default: "ALL"). Options include 'BP', 'CC', 'MF', 'ALL'.
- `--keyType`: Key type for GSEA (default: "ALIAS"). Common options include 'ENTREZID', 'SYMBOL', 'ALIAS'.
- `--nPerm`: Number of permutations for GSEA, determining the accuracy of the result (default: 10000).
- `--minGSSize`: Minimum gene set size for GSEA (default: 3).
- `--maxGSSize`: Maximum gene set size for GSEA (default: 800).
- `--pvalueCutoff`: P-value cutoff for GSEA (default: 0.05).
- `--verbose`: Verbose output (default: TRUE). Options: TRUE or FALSE.
- `--pAdjustMethod`: P-value adjustment method for GSEA (default: "none"). Options include 'BH', 'BY', 'holm', 'Bonferroni', 'none'.

## Output
The script outputs PNG files of enrichment maps that visualize significant pathways:
- `diff_corr_enrichment_map.png` for differential correlation analysis.
- `condition_1_enrichment_map.png` for Condition 1 specific analysis.
- `condition_2_enrichment_map.png` for Condition 2 specific analysis.

These images help in understanding the biological context of the gene sets involved in the study, grouped by their pathways and the strength of their enrichment scores.