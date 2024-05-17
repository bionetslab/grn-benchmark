# Methadology, Rationale and Parameter Settings

## 1. Macrophages and T-Cells
Diffcorr was performed on the Macrophages and CD_8 T cells, both of which are types of immune cells with different roles in the immune system. Macrophages are white blood cells that remove unwanted materials. They engulf and break down debris, germs, and abnormal cells. Macrophages play a vital role by detecting threats and repairing tissue damage. CD8 T cells are another type of white blood cell. They identify and eliminate virus-infected or cancerous cells. CD8 T cells release substances that destroy these harmful cells. Both cell types contribute to the immune system’s defense. The datasets contains gene expression profiles of 100 samples including 2 cell types: Macrophages and CD_8 T-cells

> Read [Downstream](https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/Downstream.md) for enrichment analysis and further information 

## 2. Golub Dataset
This dataset consist of gene expression profiles from 38 tumor samples including 2 different leukemia subtypes: 27 acute lymphoblastic leukemia (ALL) and 11 acute myeloid leukemia (AML) samples (Golub et al., 1999). 
The Affymetrix GeneChip HuGeneFL (known as HU6800) microarray contains 6800 probe-sets. 

#### Details
1. All the probe-sets with negative values in any samples were filtered out resulting in 2568 genes.
2. The expression patterns in each subtypes (ALL or AML) was then used to group the genes using the cluster.molecule function.
3. The distance mesaure (1 − correlation coefficient) was set as a cutoff value of based on the cutree function. It was set to 0.6 by the authors.
4. get.eigen.molecule and get.eigen.molecule.graph functions helped in visualising the network.
5. Finally, comp.2.cc.fdr function was used which provided the resulting pairwise differential correlations from the golub dataset. 

<img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/DiffCorr_Golub_Table.png" width="500" >

#### Results
<img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/DiffCorr_Golub_Results.png" width="900" >
Table shows the top 10 significantly differential coexpressions of all AML-specific correlations. 
For example, a correlation between D43949_at (KIAA0082) and HG4185-HT4455_at (Estrogen Sulfotransferase, Ste) was − 0.09 in ALL, and 0.98 in AML. 
The list includes genes encoding the ribosomal proteins L5, L29, L30, L37, and L37a. 
The list also contains eukaryotic translation elongation factor 1 alpha 1 (eEF1A1), which is associated with translation elongation factor activity and has oncogenic potency. 
The DiffCorr package also detetcted “switching mechanism” which are oppositely correlated pairs where, for example, 2 molecules exhibit positive correlation in one condition and negative correlation in the other condition.

## 3. AtGenExpress Leaf and Flower Samples
The dataset was downloaded using the GEOquery package. It includes microarray-based experiments measuring mRNA, genomic DNA, and protein abundance, as well as nonarray techniques such as NGS data, serial analysis of gene expression (SAGE), and mass spectrometry proteomic data.
   
A total of 34 modules in the co-expression networks with GSE5632 (flower samples) and 28 modules in the co-expression networks with GSE5630 (leaf samples) were detected.

### Results
See [Analysis](https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/Scripts/downstream_flower_leaf.R) for the code adapted to perform analysis on GSE dataset. Used the GOstats package to perform GO term enrichment analysis of the detected co-expression modules and to evaluate whether a particular molecular group is significantly over- or underrepresented. Assessed the predominant function in the biological process within the three modules

Module 1 using flower samples (GSE5632) was involved in “nucleosome assembly” within the “Biological Process” domain. Modules 2 and 3 were related to “cell proliferation” and “RNA methylation,” respectively

<img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/GO.png" width="500" >

## Parameters
- `mode="lower"` specifies that the lower triangle of the correlation matrix is used to create the graph.
- A threshold of `0.95` is used to retain edges with weights above this value in the co-expression networks.
- In hierarchical clustering, `method="pearson"` and `linkage="average"` are used.
- For `cutree`, the height parameter `h=0.4` is set to cut the hierarchical clustering tree.
- `get.eigen.molecule` performs eigenvalue decomposition (`methods="svd"`) with `n=2` components.

### Function Explanation
- Spearman correlation coefficients are calculated for correlation
- The Fruchterman-Reingold algorithm is used for network layout.
- Fast greedy community detection algorithm is applied using `fastgreedy.community`.
- Hierarchical clustering is performed with `cluster.molecule`.
- The `cutree` function is used to cut the hierarchical clustering tree at a specified height.
- Eigen decomposition is performed on the correlation matrices for community detection using `get.eigen.molecule`.
- Visualization graphs are plotted using `plot`.
- Module assignments for each gene are written to files using `write.modules`.
- Differential correlation analysis is performed using the `plotDiffCorrGroup` and `comp.2.cc.fdr` functions from the `DiffCorr` package.
- The `plotDiffCorrGroup` function is used to visualize differential correlations between groups.
<p align="right">(<a href="#readme-top">back to top</a>)</p>
