# Enrichment Analysis
The DiffCorr package offers a list of pairwise differentially correlated molecules. the file also includes r-values, p-vales, and fdr values. By identifying the pair corresponding to the maximum absolute difference of r-values (r1-r2), indicating the most differentially correlated pair, enrichment analysis can be perfomed to obtain further useful information

## gProfiler
gProfiler (Gene List Profiling) is a web server and tool suite that is used for functional annotation and enrichment analysis of gene lists. It takes sets of genes and maps them onto known biological pathways, ontologies and other annotation categories. Enrichment analysis is carried out to identify over-represented biological terms as well as pathways in a given gene list. Besides, it annotates genes with information obtained from various databases such as Gene Ontology (GO), KEGG, Reactome, WikiPathways among others.

### Step 1: Sort the differentially correlated pairs of genes in descending order of (r1 - r2) values
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/rSortedList.png" width="700" >
</div>

### Step 2: Upload the list of most differentially correlated pairs of genes to gProfiler Web
Visit  [gProfiler Web App](http://biit.cs.ut.ee/gprofiler/) and upload required list of genes for analysis
For the purpose of explaining the analysis better, only the first pair of genes were uploaded. _More pair can be uploaded for extensive analysis._
After the query has been filled out, Press _Run query_ button
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/uploadGenes.png" width="700" >
</div>

### Step 3: Visulaze different ascpects of analysis as per requirements 
##### Overview
The first tab shows the significance of different functional categories across various databases. -log₁₀(padj) signifies the enrichment: higher the log value, more significant enrichment.

The table below explains detailed information about the enriched terms. It includes __Source__ displaying the database or ontology from which the term is sourced from, __Term ID__ which is a unique identifier for the term within its source database, __term name__ and __adjusted p-value__.

<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/overview.png" width="700" >
</div>

##### Detailed Results
Provides valuable insights into the functional roles of molecules in cellular processes. It would be better if there could be additional experiments carried out such as protein-protein interaction assays and functional tests so as to establish more about the functions that these genes play in various pathways during different biological processes. The results shown here are generated from diverse biological resources and databases which are popularly employed in genomic science concerning functional annotation as well as pathway analysis.

For instance, the Gene Ontology Molecular Function [GO:MF] highlights how RBM3 and C1QBP genes play important roles in cells. They are good at binding to ribosomes and ribonucleoprotein complexes, helping to make proteins and process genetic information. They're also involved in binding to kininogen, a protein that helps regulate blood clotting and control inflammation in the body.

<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/detailedResults.png" width="700" >
</div>

##### GO Context
gProfiler also provides a context network containing nodes with directed edges that reflects the ontology structure of biological processes. It is more like a visual summary of the enriched biological processes, helping to identify major functional themes and their relationships.
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/detailedResults.png" width="700" >
</div>

Positive Regulation of Gene Expression (GO:0010628): This term suggests that RBM3 and C1QBP are involved in processes that increase the expression of specific genes. This could be through various mechanisms such as enhancing transcription, increasing mRNA stability, or promoting translation.
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/gopb_1.png" width="700" >
</div>

Regulation of Nitrogen Compound Metabolic Process (GO:0051171): This term indicates that these genes play roles in the regulation of metabolic processes involving nitrogen-containing compounds. This includes amino acids, nucleotides, and other nitrogenous substances.
<div align="center">
    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/gopb_2.png" width="700" >
</div>

> Note: More pairs of genes can be analysed together
> <div align="center">
>    <img src="https://github.com/aparnaullas97/grn-benchmark/blob/main/src/diffcorr/ImageResouces/cumulative.png" width="700" >
> </div>

