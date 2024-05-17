# Downstream Analysis
The downstream analysis was performed on three datasets: leaf and flower samples from the AtGenExpress development (Accession: GSE5630 and GSE5632, respectively), the Golub dataset and the Macrophages dataset. Here we go through the details of downstream analysis performed on the Arabidopsis dataset

### Heatmaps
After correlation of both the datasets, spatstat package was used to plot the heatmaps
Horizontal and vertical of the heatmaps of the gene expression correlation matrices show the probe set identifiers in each experiment.
Pink = positive correlation, blue = negative correlation between the two probe
sets.
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/heatmap.png" width="400" >
</div>

### Construction of the co-expression network
Next, using the igraph package, the co-expression graph was constructed and visualised. The threshold value, rs was set as ≥ 0.95. 
Here, nodes are the probe sets and the edges mean that the correlation coefficients are over 0.95 between the connected nodes. The nodes with degree > 20 were colored magenta and those that below, green.
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/coexpress.png" width="400" >
</div>
In order to visualize better, the graphs were exported as 'gml' files which help networks become more interactive. These files then can be imported into Cytoscape for further evaluation

### Graph Clustering
In order to find gene co-expression modules, igraph package implementation named Fast Greedy graph clustering algorithm was used. This algorithm runs in linear time O(n log2n), on a network with n vertices and reduces computation time.
After which each module member could be accessed using the membership functionality.

### Gene Ontology Enrichment Analysis
Finally to access cluster fidelity, Gene Ontology (GO) Biological Process (BP) enrichment analysis was performed. In this analysis, a group of genes undergoes testing to determine whether there is an over-representation of particular biological processes as defined by Gene Ontology terms.

Breaking down each column:
- **GOBPID**: This is the identification code assigned to each Gene Ontology Biological Process term.
- **Pvalue**: This value is the statistical significance level, showing how likely it is to observe the given number of genes associated with the GO term by random chance.
- **OddsRatio**: This quantifies the likelihood of observing the number of genes linked to the GO term in the input gene set compared to what would be expected by chance.
- **ExpCount**: This helps understand the expected count of genes associated with the GO term based on the background distribution.
- **Count**: This denotes the actual count of genes identified from the input gene set that are associated with the specific GO term.
- **Size**: This shows the total number of genes annotated to the GO term in the reference database.
- **Term**: This is the name or description of the GO BP term, describing the biological process it represents.

A low p-value indicates that it's highly unlikely the connection between the genes we're examining and the GO term occurred randomly. Conversely, a high odds ratio suggests a robust connection between the genes and the GO term.

For instance, in the first row, when we look at the GO BP term "ribosome biogenesis," its p-value of 0.000 signals a significant enrichment. The odds ratio of 39.746 emphasizes a strong link between the genes we're studying and this biological process. Specifically, we find 9 genes associated with "ribosome biogenesis" in our set, whereas we'd expect 0 based on the background distribution. Moreover, there are a total of 316 genes annotated to "ribosome biogenesis" in the reference database.
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/GOEA.png" width="400" >
</div>
