
## Execution Commands for Reference Datasets

This section provides detailed commands for executing the DGCA tool with reference datasets both using Docker and directly via Rscript.

### With Docker

Execute the DGCA analysis using Docker for containerized environment management:

**500 Genes:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data dgca-tool dgca.R --input_file_1 /data/500_genes/out_CD8_exhausted.tsv --input_file_2 /data/500_genes/out_Macrophages.tsv --output_path /data/500_genes
```

**1000 Genes:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data dgca-tool dgca.R --input_file_1 /data/1000_genes/out_CD8_exhausted.tsv --input_file_2 /data/1000_genes/out_Macrophages.tsv --output_path /data/1000_genes
```

**2500 Genes:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data dgca-tool dgca.R --input_file_1 /data/2500_genes/out_CD8_exhausted.tsv --input_file_2 /data/2500_genes/out_Macrophages.tsv --output_path /data/2500_genes
```

**Full Input:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data dgca-tool dgca.R --input_file_1 /data/full_input/out_CD8_exhausted.tsv --input_file_2 /data/full_input/out_Macrophages.tsv --output_path /data/full_input
```

### Without Docker

Execute the DGCA analysis using Rscript directly for various gene datasets as follows:

**500 Genes:**
```bash
Rscript dgca.R --input_file_1 ../../reference_datasets/500_genes/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/500_genes/out_Macrophages.tsv --output_path ../../reference_datasets/500_genes
```

**1000 Genes:**
```bash
Rscript dgca.R --input_file_1 ../../reference_datasets/1000_genes/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/1000_genes/out_Macrophages.tsv --output_path ../../reference_datasets/1000_genes
```

**2500 Genes:**
```bash
Rscript dgca.R --input_file_1 ../../reference_datasets/2500_genes/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/2500_genes/out_Macrophages.tsv --output_path ../../reference_datasets/2500_genes
```

**Full Input:**
```bash
Rscript dgca.R --input_file_1 ../../reference_datasets/full_input/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/full_input/out_Macrophages.tsv --output_path ../../reference_datasets/full_input
```
