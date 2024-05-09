## SCANet:

### Methodology:
1. Representative cells are created to reduce dimensionality and computational complexity by grouping neighboring cells in the preprocessed data, which is an expression matrix with rows representing individual cells and columns representing genes
2. It utilizes modified Weighted Gene Co-expression Network Analysis (WGCNA) to infer modules of strongly co-expressed genes.
3. SCANet allows users to explore and visualize gene co-expression networks (GCNs), helping them identify specific modules of interest for further investigation.
4. Gene regulatory networks (GRN) of the modules of interest are inferred using GRNBoost2 algorithm. RcisTarget verifies accuracy through cis-regulatory TF-binding motif enrichment analysis.
5. By comparing the gene networks with drug databases, SCANet suggests existing drugs that could potentially influence these gene networks.

### Rationale:
SCANet is designed to simplify scRNA-seq data analysis by combining established bioinformatics techniques and algorithms in a comprehensive platform. It streamlines the process of constructing co-expression and regulatory networks while also linking these networks to potential drug repurposing opportunities, helping researchers uncover critical regulatory mechanisms and therapeutic targets in different disease contexts.




### Parameters:
- Number of representative cells per cluster.
- correlation_cutoff, which Specifies the minimum threshold for co-expression strength that defines whether an edge should be included in the network. Edges with co-expression values below this threshold will be excluded from the network.
