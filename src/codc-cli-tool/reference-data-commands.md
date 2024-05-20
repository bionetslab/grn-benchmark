
## Execution Commands for Reference Datasets

### With Docker

**500 Genes:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data codc-tool codc --input_file_1 /data/500_genes/out_CD8_exhausted.tsv --input_file_2 /data/500_genes/out_Macrophages.tsv --output_path /data/500_genes
```

**1000 Genes:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data codc-tool codc --input_file_1 /data/1000_genes/out_CD8_exhausted.tsv --input_file_2 /data/1000_genes/out_Macrophages.tsv --output_path /data/1000_genes
```

**2500 Genes:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data codc-tool codc --input_file_1 /data/2500_genes/out_CD8_exhausted.tsv --input_file_2 /data/2500_genes/out_Macrophages.tsv --output_path /data/2500_genes
```

**Full Input:**
```bash
docker run --rm -v $(pwd)/../../reference_datasets:/data codc-tool codc --input_file_1 /data/full_input/out_CD8_exhausted.tsv --input_file_2 /data/full_input/out_Macrophages.tsv --output_path /data/full_input
```

### Without Docker

**500 Genes:**
```bash
pdm run cli codc --input_file_1 ../../reference_datasets/500_genes/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/500_genes/out_Macrophages.tsv --output_path ../../reference_datasets/500_genes
```

**1000 Genes:**
```bash
pdm run cli codc --input_file_1 ../../reference_datasets/1000_genes/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/1000_genes/out_Macrophages.tsv --output_path ../../reference_datasets/1000_genes
```

**2500 Genes:**
```bash
pdm run cli codc --input_file_1 ../../reference_datasets/2500_genes/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/2500_genes/out_Macrophages.tsv --output_path ../../reference_datasets/2500_genes
```

**Full Input:**
```bash
pdm run cli codc --input_file_1 ../../reference_datasets/full_input/out_CD8_exhausted.tsv --input_file_2 ../../reference_datasets/full_input/out_Macrophages.tsv --output_path ../../reference_datasets/full_input
```
