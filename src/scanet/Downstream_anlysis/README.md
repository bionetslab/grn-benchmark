
# Downstream Analysis

This script conducts downstream analysis on the inferred gene regulatory network, utilizing GO enrichment analysis and KEGG pathway enrichment analysis.

Initially, the network is divided into two distinct networks based on the respective conditions: one for the exhausted condition and the other for the macrophage condition. Subsequently, Go enrichment analysis and KEGG analysis are conducted on each of these networks.


## Requirements

- Python 3.6+
- argparse
- gprofiler
- gseapy
- pandas
- os

Ensure you have Python 3.6+ and install the required packages:

## Input

- `Path to the network.tsv` file, which is a tab-separated file that contains all edges with the columns target, regulator, condition and weight 

- `output path` : String indicating the path where the output will be stored. The directory doesn't have to exist beforehand; it will be automatically created if it doesn't


## Output

The output files will be saved in the specified directory. These will include:

- **Go_enriched_terms.tsv:** The results of the GO enrichment analysis
  
- **Top_go_enriched_terms.tsv:** Here, the top enriched GO terms for each condition, based on the ratio between the `intersection size` and `term size`.
  
- **KEGG_enriched_pathways.tsv** The results of the KEGG pathway enrichment analysis.
  
- **Top_KEGG_enriched_pathways.tsv:** The top enriched KEGG pathways, based on the `Combined Score`.


## Usage

Provide an input network file in TSV format and specify an output directory for results:

```bash
python downstream_analysis.py --c <input_network.tsv> --o <output_directory>
```

### Arguments

- `--c`: Path to the input network file (TSV format).
- `--o`: Directory where output files will be saved.

## Example

```bash
python downstream_analysis.py --c "../../network.tsv" --o "../../output"
```
