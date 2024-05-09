# SCANeT on the reference datasets

This Python script offers features for inferring gene co-expression networks (GCNs) and gene regulatory networks (GRNs), as well as identifying potential drug candidates using the provided reference datasets. Initially, the script reads the data, converts it to the h5ad format, and visualizes it. Subsequently, representative samples are obtained by averaging neighboring samples. Based on the correlation between genes, modules of genes are generated. The most two intriguing modules, exhibiting the highest correlation with the two conditions, are selected. The GRNs for these selected modules are inferred. Finally, potential drug candidates are identified.


## Requirements

- Python 3.6+
- scanet
- anndata
- scanpy
- pandas
- os

## Method parameters

- `N` : Number of representative cells/samples per cluster.
- `correlation_cutoff` : Specifies the minimum threshold for co-expression strength that defines whether an edge should be included in the network. Edges with co-expression values below this threshold will be excluded from the network.

## Usage

The script can be executed from the command line with appropriate arguments:

```bash
python scanet_on_ref_data.py --c <input_file_exhausted> --d <input_file_macrophages> --o <output_directory> --n <num_rep_cells> --t <correlation_cutoff>
```


## Inputs

- `input_file_exhausted`: Path to the first input file which is the exhausted.tsv file.
- `input_file_macrophages`: Path to the second input file which is the Macrophages.tsv file.
- `output_directory`: Directory where the output will be stored.
- `num_rep_cells` : Number of representative cells. It has the default value 75.
- `correlation_cutoff` : Minimum threshold for co-expression strength. It has the default value 0.85.

## Outputs

In the provided directory, you will find the following :

- `network_exhausted` : A tsv file which represents the GRN of the module, which has the highest correlation to exhausted condition
- `network_macrophages` : A tsv file which represents the GRN of the module, which has the highest correlation to macrophages condition
- `outs_scan` : A folder which has the plots of gcn and grn with/without drug candidates for the two conditions. This is not required, but I added them. 

## Example

To run the script with example data:

```bash
python scanet_on_ref_data.py --c "data\out_CD8_exhausted.tsv" --d "data\out_Macrophages.tsv" --o "output/" --n 75 --t 0.85
```

Replace `data\out_CD8_exhausted.tsv`, `data\out_Macrophages.tsv`, and `output/` with your actual file paths and desired output location.

The last two arguments `--num_rep_cells` and `--correlation_cutoff` are not required.
