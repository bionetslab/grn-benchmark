# Downstream analysis
A suitable downstream analysis of the network is to analyze for Gene Ontology (GO) terms enrichment and KEGG pathways enrichment of the network genes.

The list of genes for every condition in the network is analyzed for GO terms enrichment and KEGG pathways enrichment using `clusterProfiler` package (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).

### Dependencies:
- R
- BiocManager
- clusterProfiler
- org.Hs.eg.db
- tidyverse
- data.table
- optparse
- Rscript


### Examplery Execution Instructions:
```
Rscript downstream.R -n network.tsv -o ./
```

### Input file format specification:
`--network` or `-n`: Path to tab-separated file that contains the fowllowing columns:
- `target`: Targeted genes
- `regulator`: Regulatory genes
- `condition`: Condition that the edge belongs to
- `weight`: Weight of the edge

### Output file format specification:
Two tab-separated output files stored at `--output_path`:
- `GO_enriched_terms.tsv` containing the results of the GO terms enrichment analysis. The file contains the following columns: 
   - `ID`: GO term ID
   - `Description`: Description of the GO term
   - `GeneRatio`: Ratio of the number of genes in the gene list that are associated with the GO term to the total number of genes associated with the GO term
   - `BgRatio`: Ratio of the number of genes in the background gene list that are associated with the GO term to the total number of genes associated with the GO term
   - `pvalue`: P-value of the enrichment analysis
   - `p.adjust`: Adjusted p-value of the enrichment analysis
   - `qvalue`: Q-value of the enrichment analysis
   - `geneID`: Gene IDs associated with the GO term
   - `Count`: Number of genes associated with the GO term
   - `condition`: Condition that the GO term belongs to
  
- `KEGG_enriched_pathways.tsv` containing the results of the KEGG 
pathways enrichment analysis. The file contains the following columns: 
    - `category`: Category of the pathway
    - `subcategory`: Subcategory of the pathway
    - `ID`: Pathway ID
    - `Description`: Description of the pathway
    - `GeneRatio`: Ratio of the number of genes in the gene list that are associated with the pathway to the total number of genes associated with the pathway
    - `BgRatio`: Ratio of the number of genes in the background gene list that are associated with the pathway to the total number of genes associated with the pathway
    - `pvalue`: P-value of the enrichment analysis
    - `p.adjust`: Adjusted p-value of the enrichment analysis
    - `qvalue`: Q-value of the enrichment analysis
    - `geneID`: Gene IDs associated with the pathway
    - `Count`: Number of genes associated with the pathway
    - `condition`: Condition that the pathway belongs to