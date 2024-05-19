The study outlines a systematic approach for conducting Differential Network Analysis (DiNA) to identify differentially co-expressed modules (DCMs) in gene expression data. Here’s a summarized view of the methodology, presented as a series of steps:

1. **Identification of Gene Modules**:
   Gene modules are identified as the initial step. Each module must contain at least three genes and can be sourced from pre-existing databases or determined through clustering techniques. Notably, modules may overlap, allowing genes to appear in more than one module.

2. **Construction of Gene Expression Matrices**:
   For each identified module, a gene expression matrix \(X(gm)\) is constructed for various groups, such as healthy versus diseased individuals. The expression levels in these matrices may be represented as integer counts from microarray data or continuous values from sequencing technologies.

3. **Development of Similarity Matrices**:
   A similarity matrix \(S(gm)\) is created for each group to analyze the interactions within each module. The matrix is symmetric and its size depends on the number of genes in the module. The matrix elements are calculated using different correlation measures like Pearson, Spearman, or Kendall, or other methods like partial correlation or mutual information.

4. **Hypothesis Testing on Network Structures**:
   The core of the analysis involves hypothesis testing to compare network structures across groups, represented by the similarity matrices. The goal is to determine if there are significant differences in the module interactions between groups, comparing matrices \(S(g1m)\) and \(S(g2m)\).

5. **Implementation of the p-Norm Difference Test (PND)**:
   Among the statistical tests utilized, the p-norm difference test (PND) is introduced as a new approach. This test is designed to effectively identify significant differences in network structures while maintaining a low rate of false positives.

6. **P-Value Calculation via Non-Parametric Permutation**:
   P-values are calculated using a non-parametric permutation method, ideal for data with complex dependencies, as it accurately reflects the structured nature of similarity matrices.

7. **Software Implementation**:
   The entire methodology and associated tests are implemented in the R package `discoMod`, which provides tools for clustering gene modules, testing for differential co-expression, and visualizing results. This package streamlines the process for researchers, facilitating easier analysis and interpretation of gene expression data.
