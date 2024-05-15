# chNet_BioNet

### Differential network analysis by simultaneously considering changes in gene interactions and gene expression
#### Jia-Juan Tu, Le Ou-Yang, Yuan Zhu, Hong Yan, Hong Qin, Xiao-Fei Zhang
<br>

```chNet_BioNet``` is an advanced wrapper tool to differential network analysis by simultaneously analyzing differential gene interaction and differential gene expression. This robust tool identifies key biological differences across various conditions, like ```Macrophages``` and ```CD8_exhausted```. Here, I tried to reimplement the [chNet](https://github.com/Zhangxf-ccnu/chNet) and also apply futher down stream analysis upon it.

## Installation 
Before delving into the code, please check these prerequisites.
- [ ] Install RStudio 
- [ ] Install R 4.3.3 and define the PATHs
- [ ] Install Docker

*Note:* I just tested my commands in the `Windows` environment. I haven't tested the feasibility of the commands in `Linux` or `OS X`.

<br>

To run the project, you have two solutions. 
- **Docker**
- **Rscript Command Line**
<br>

## Project's Structure 
```md
chNet_BioNet
├── .git
│   └── (git files and directories)
│   │    
├── docker
│   ├── Dockerfile
│   ├── downstream.R
│   ├── input_output.R
│   ├── replicate.R
│   ├── output
│   │   ├── downstream
│   │   │   ├── 500_genes
│   │   │   │   └── (PDF files)
│   │   │   ├── 1000_genes
│   │   │   │   └── (PDF files)
│   │   │   ├── 2500_genes
│   │   │   │   └── (PDF files)
│   │   │   ├── GSE13159.AML
│   │   │   │   └── (PDF files)
│   │   │   └── TCGA.BRCA
│   │   │       └── (PDF files)
│   │   ├── input_output
│   │   │   ├── 500_genes
│   │   │   │   └── (PDF files)
│   │   │   ├── 1000_genes
│   │   │   │   └── (PDF files)
│   │   │   ├── 2500_genes
│   │   │   │   └── (PDF files)
│   │   │   ├── GSE13159.AML
│   │   │   │   └── (PDF files)
│   │   │   └── TCGA.BRCA
│   │   │       └── (PDF files)
│   │   └── replicate
│   │       ├── 500_genes
│   │       │   └── (PDF files)
│   │       ├── 1000_genes
│   │       │   └── (PDF files)
│   │       ├── 2500_genes
│   │       │   └── (PDF files)
│   │       ├── GSE13159.AML
│   │       │   └── (PDF files)
│   │       └── TCGA.BRCA
│   │           └── (PDF files)
│   ├── reference_datasets
│   │   ├── 500_genes
│   │   │   ├── out_CD8_exhausted.tsv
│   │   │   └── out_Macrophages.tsv
│   │   ├── 1000_genes
│   │   │   ├── out_CD8_exhausted.tsv
│   │   │   └── out_Macrophages.tsv
│   │   └── 2500_genes
│   │       ├── out_CD8_exhausted.tsv
│   │       └── out_Macrophages.tsv
├── R
│   └── (same files and directories' as docker directory.)
│   │    
└── README.md
```

The ```chNet_BioNet``` project includes a docker directory that contains R scripts for network analysis and their corresponding output directories, each holding PDF results for various datasets. The reference_datasets folder houses essential data files. The R directory mirrors the docker setup for consistency. This structure ensures easy navigation and use of the tools within Docker or R environments.

### Docker
At first, by using this command, you can build the dockerfile, and continue to run each part of the project with these docker commands.
```Dockerfile
docker build -t r-scripts .
```
Having generated the dockerfile, you can see now the general command for running my project. As mentioned in the **tree**, I have three types of output in my project.
- Input/Output Analysis
- Downstream Analysis
- Replication Analysis
<br>

#### General command line
```Dockerfile
docker run \
  -v "${PWD}/output:/usr/src/app/output" \
  -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" \ # Optional: only for external datasets
  --type of script r-scripts Rscript YourScript.R \
  --type <internal|external> \
  --dataset "DATASET_NAME" \ # Optional: used for internal datasets
  --input1 "FILE_PATH_1" \ # Optional: used for the first input file in external datasets
  --input2 "FILE_PATH_2" \ # Optional: used for the second input file in external datasets
  --output "OUTPUT_DIRECTORY" \
  --resultFile "RESULT_FILE_NAME"
```

##### Parameter Explanations
I have two important types of dataset. The first one is ```internal``` which is used for chNet package's dataset, including **TCGA.BRCA** and **GSE13159.AML** datasets. And the second type of dataset is ```external```, which shows the ```reference dataset```. Here all of the arguments of the docker command is explained with details.

- *--v "${PWD}/output:/usr/src/app/output"*: Maps the local output directory to the corresponding directory inside the Docker container.

- *--v "${PWD}/reference_datasets:/usr/src/app/reference_datasets*: This parameter is only used for external datasets. It maps the reference_datasets directory on the host to the container.

- *--type of script*: You can choose between **input_output.R**, **downstream.R** and **replicate.R** based on desired output.

- *--type*: This parameter indicates whether the dataset is **internal** or **external**.

- *--dataset*: Specifies the name of the dataset when using **internal** datasets. This parameter is used to select which built-in dataset to process.

- *--input1* and *--input2*: These parameters specify the file paths for the first and second input files respectively when using **external** datasets.

- *--output*: Specifies the output directory path inside the container where result files will be stored.

- *--resultFile*: Defines the name of the output file to be generated.

Now, let's dive into running the code. 

<br>

### Minimal Working Examples

#### Input/Output Analysis
In this part, by using each of these commands, you will be able to generate your desired output for your prefered dataset. The output file is in the form of ```type_network.tsv``` which includes ```regulator```, ```target```, ```condition``` and ```weight```.

```Dockerfile
docker run -v "${PWD}/output:/usr/src/app/output" r-scripts Rscript input_output.R internal "TCGA.BRCA" "output/input_output" "TCGA_BRCA_network.tsv"
docker run -v "${PWD}/output:/usr/src/app/output" r-scripts Rscript input_output.R internal "GSE13159.AML" "output/input_output" "GSE13159.AML_network.tsv"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript input_output.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "output/input_output" "gene500_network.tsv"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript input_output.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "output/input_output" "gene1000_network.tsv"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript input_output.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "output/input_output" "gene2500_network.tsv"
```

#### Reproduction Analysis
For reproduction of the paper's figures, you can run each of these docker commands accordingly. For each of these dataset, you can find the output directory which has several pdf files for the specified dataset.

```Dockerfile
docker run -v "${PWD}/output:/usr/src/app/output" r-scripts Rscript replicate.R internal "TCGA.BRCA" "Basal" "LumA" "output/replicate/TCGA.BRCA" "TCGA.BRCA"
docker run -v "${PWD}/output:/usr/src/app/output" r-scripts Rscript replicate.R internal "GSE13159.AML" "cancer" "normal" "output/replicate/GSE13159.AML" "GSE13159.AML"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript replicate.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/500_genes" "gene500"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript replicate.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/1000_genes" "gene1000"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript replicate.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/2500_genes" "gene2500"
```

#### Downstream Analysis
For further analysis upon the paper, several figures are also provided. By running each of these, the same logic with the previous parts, you will be able to have pdf files as output.

```Dockerfile
docker run -v "${PWD}/output:/usr/src/app/output" r-scripts Rscript downstream.R internal "TCGA.BRCA" "Basal" "LumA" "output/downstream/TCGA.BRCA" "TCGA.BRCA"
docker run -v "${PWD}/output:/usr/src/app/output" r-scripts Rscript downstream.R internal "GSE13159.AML" "cancer" "normal" "output/downstream/GSE13159.AML" "GSE13159.AML"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript downstream.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/500_genes" "gene500"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript downstream.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/1000_genes" "gene1000"
docker run -v "${PWD}/output:/usr/src/app/output" -v "${PWD}/reference_datasets:/usr/src/app/reference_datasets" r-scripts Rscript downstream.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/2500_genes" "gene2500"
```
<br>

### Rscript Command Line
According to the [chNet](https://github.com/Zhangxf-ccnu/chNet) R package on GitHub, it should be installed as below.
```R
# Step 1. Install the devtools package. Invoke R and then type
install.packages("devtools") 

# Step 2. Load the devtools package.
library("devtools") 

# Step 3. Install the chNet package from GitHub.
install_github("Zhangxf-ccnu/chNet", subdir="pkg")
```

#### Running Commands in RStudio's Terminal on Windows:
Go to the "Tools" menu and select ```Terminal``` to open the Terminal tab. Navigate to the directory where you want to run the R commands (if needed) using the cd command.
Make sure ```R (version 4.3.3 in your case)``` and RStudio are installed on your Windows system before executing these commands.

Having all the libraries installed in your system, now you are enabled to run these R-commands. Just be aware of the fact, all of the parameters here are as same as the docker commands. So, for further explanation of each parameter, you can take a look again to docker part's command's explanation.

<br>

### Minimal Working Examples
#### Input/Output Analysis
```R
Rscript input_output.R internal "TCGA.BRCA" "output/input_output" "TCGA_BRCA_network.tsv"
Rscript input_output.R internal "GSE13159.AML" "output/input_output" "GSE13159.AML_network.tsv"
Rscript input_output.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "output/input_output" "gene500_network.tsv"
Rscript input_output.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "output/input_output" "gene1000_network.tsv"
Rscript input_output.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "output/input_output" "gene2500_network.tsv"
```

#### Reproduction Analysis
```R
Rscript replicate.R internal "TCGA.BRCA" "Basal" "LumA" "Output/replicate/TCGA.BRCA" TCGA.BRCA
Rscript replicate.R internal "GSE13159.AML" "cancer" "normal" "output/replicate/GSE13159.AML" GSE13159.AML
Rscript replicate.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/500_genes" gene500
Rscript replicate.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/1000_genes" gene1000
Rscript replicate.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/replicate/2500_genes" gene2500
```

#### Downstream Analysis
```R
Rscript downstream.R internal "TCGA.BRCA" "Basal" "LumA" "output/downstream/TCGA.BRCA" TCGA.BRCA
Rscript downstream.R internal "GSE13159.AML" "cancer" "normal" "output/downstream/GSE13159.AML" GSE13159.AML
Rscript downstream.R external "reference_datasets/500_genes/out_CD8_exhausted.tsv" "reference_datasets/500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/500_genes" gene500
Rscript downstream.R external "reference_datasets/1000_genes/out_CD8_exhausted.tsv" "reference_datasets/1000_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/1000_genes" gene1000
Rscript downstream.R external "reference_datasets/2500_genes/out_CD8_exhausted.tsv" "reference_datasets/2500_genes/out_Macrophages.tsv" "CD8 Exhausted T-cells" "Macrophages" "output/downstream/2500_genes" gene2500
```
<br>

## Output
- ### Input/Output Analysis
* `network.tsv`: tab-separated file that contains all edges (row-wise) with the following columns:
     * First column `target`: Target of the edge
     * Second column `regulator`: Source of the edge
     * Third column `condition`: Whether `Differential` or `Non-differential`.
     * Fourth column `weight`: Weight of the edge, in my case all is `1`, because the weighted chNet has an [error](#error).

 *Note:* In this [part](#Agorithm), I briefly explained how to find Differential nodes and interactions.
 
<br>

## Method and Result
This study introduces an advanced differential network analysis method designed to identify changes in gene interactions alongside in gene expression levels. This approach, highlighted in the paper ```Differential network analysis by simultaneously considering changes in gene interactions and gene expression``` integrates hypothesis tests to evaluate changes in partial correlations between gene pairs with changes in their expression. The methodology employs an optimization framework to enforce hierarchical constraints which improves the interpretability and accuracy of the results. The proposed method has been validated through simulations and applied to case studies including breast cancer subtypes and acute myeloid leukemia, demonstrating its superiority over existing methods. The algorithm can be seen here.

<br>

<a name="Agorithm"></a>
### General Agorithm
#### Input: Gene Expression Datasets: X(1) and X(2)
- **Tuning Parameter:** (`lambda = 2.825`): controls the level of sparsity and a threshold parameter.

#### Expected Output:
- **Differential Edges**
- **Differentially Expressed Genes**

#### Procedure:
- **Compute Test Statistics** for partial correlations (`t_ij`) (quantifies the loss) and expression levels (`z_i`) (impose the sparsity).
- **Calculate Thresholds** (`lambda_i` and `lambda_ij`) to identify significant changes.
- **Determine Significant Interactions** and **gene expressions** using the calculated thresholds.

#### Results:
- `Genes` and `Edges` that are bigger than the defined threshold, (`lambda_i` and `lambda_ij`) respectively. In this case `Diff`, otherwise `Non-diff`.
  
<br>

### Figure 1: Reproduction Analysis

<table>
  <tr>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/Difference_of_Partial_Correlation_gene500.png" alt="Difference of Partial Correlation" width="300"/></td>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/Differential_network_graph_gene500.png" alt="Differential Network Graph" width="300"/></td>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/Hub_genes_boxplot_gene500.png" alt="Hub Genes Boxplot" width="300"/></td>
  </tr>
  <tr>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/Expression_data_CD8_exhausted_gene500.png" alt="Expression Data CD8 Exhausted" width="300"/></td>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/Expression_data_Macrophages_gene500.png" alt="Expression Data Macrophages" width="300"/></td>
  </tr>
</table>


This comprehensive visualization (Figure 1) captures the essence of differential network analysis by showcasing several key aspects as below:

- **1. Difference of Partial Correlation**: The top-left plot visualizes the difference in partial correlations between gene pairs under two conditions. Each dot represents a significant change in correlation, aiding in identifying potential areas of interest for further analysis.

- **2. Differential Network Graph**: Displayed in the top-center, this network graph illustrates the type of differentially expressed genes to be whether ```Diff Gene``` and ```Non-diff Gene```, differentiating with **Red** and **Purple** respectively. Nodes (genes) are sized according to their degree of connectivity, with hub genes highlighted.

- **3. Hub Gene Boxplot**: The top-right plot provides a boxplot of gene expression levels for hub genes identified in the differential network. Different colors represent different conditions (e.g., CD8 exhausted and Macrophages), illustrating variations in expression that correlate with gene significance in the network analysis as well as the scatters, showing with dots.

- **4. Expression Data Heatmaps**: The two large plots at the bottom show heatmaps of gene expression data for the each conditions studied. These heatmaps provide a visual representation of expression levels across different samples, categorized by expression intensity from low (green) to high (red), indicating comparison of expression patterns between these two conditions (e.g., CD8 exhausted and Macrophages).

<br>

### Figure 2: Down Stream Analysis
<table>
  <tr>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/network_clustering_gene500.png" alt="Network Clustering" width="300"/></td>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/network_histogram_gene500.png" alt="Network Histogram" width="300"/></td>
    <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/hub_gene_correlation_gene500.png" alt="Hub Gene Correlation" width="300"/></td>
  </tr>
  <tr>
  <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/hub_gene_boxplot_gene500.png" alt="Hub Gene Expression Boxplot" width="300"/></td>
  <td><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/hub_gene_heatmap_gene500.png" alt="Hub Gene Heatmap" width="300"/></td>
  </tr>
</table>

This figure presents a comprehensive downstream analysis of gene interactions and expressions using the [chNet](https://github.com/Zhangxf-ccnu/chNet) package:

- **1. Clustering of Differential Network Graph**:
The differential network's clustering according to gene interactions is shown in the top-left graph. Every node represents a gene, with colors determined by the community. This clustering aids in locating connected gene groupings that might have comparable functional characteristics.

- **2. Degree Distribution of Network Nodes**:
The top-center histogram visualizes the degree distribution within the differential network. This plot shows the frequency of genes having a specific number of connections, which is crucial for identifying highly connected hub genes within the network.

- **3. Hub Gene Expression Correlation Plot**:
The top-right plot uses a correlation matrix to show the relationship between hub genes. Each circle’s size and color intensity represent the strength and direction of the correlation, providing insights into how hub genes may co-regulate each other.

- **4. Hub Gene Expression Boxplot**:
The bottom-left boxplot displays expression levels of hub genes across two conditions, **CD8 exhausted** and **Macrophages**. This visualization helps compare gene expression variability and central tendency between groups, offering insights into gene behavior under different biological states.

- **5. Heatmap of Gene Expression**:
The bottom-right heatmap illustrates the **expression levels** of significant genes across samples, grouped by their biological condition. The colors indicate **z-score normalized expression levels**, with hierarchical clustering applied to both rows (genes) and columns (samples). This heatmap effectively shows patterns of expression across conditions and identifies genes with similar expression profiles.

<br>

### Figure 3: Gene Antology Enrichment
<table>
 <tr>
      <td colspan="2"><img src="https://github.com/sinarazi/chNet_BioNet/blob/main/documents/pics/gene_ontology_enrichment_gene500.png" alt="Gene Ontology Enrichment" width="600"/></td>
  </tr>
</table>

This bar chart illustrates the results of the Gene Ontology (GO) enrichment analysis, using ```clusterProfiler```, performed on the top hub genes identified from the differential network analysis. The analysis focuses on biological processes, highlighting those that are significantly enriched among these hub genes. Each bar represents a specific **biological process**, with the height indicating the **number of connected genes** involved in that process. This provides insights into the functional roles that these hub genes may play in cellular processes.

<br>

## Discussion 
### Differential Network Analysis
This approach, ```chNet```, is presented by ```Tu et al```, and uses hierarchical constraints to ensure that only significant, biologically relevant gene interactions are considered. By successfully filtering out indirect effects, this method guarantees that modifications in gene interactions are supported by consistent differences in gene expression level. Because this technique only takes into account interactions supported by observable changes in expression, it exhibits better accuracy and biological relevance. 
 
### Downstream Analysis
```chNet_BioNet``` adds powerful downstream analysis tools and enhanced visualization capabilities, going beyond the initial capabilities. Deeper understanding and hypothesis testing in clinical studies are made possible by these improvements, which enable a more thorough analysis and display of gene interactions and expression patterns under different conditions.
<br>

<a name="error"></a>
#### Package's Error
The ```chNet paper``` has two algorithms. The first one is the normal chNet algorithm which has the weight 1. The second algorithm uses sub-sampling to generate **weighted differential network**. Nevertheless, when it comes to using the second algorithm in the R package, I faced this issue, showing in below.

```R
> library(chNet)
> data("TCGA.BRCA")
> result = chNet(TCGA.BRCA$X,TCGA.BRCA$group,subsampling = TRUE, R= 20,
+                lambar = 2.825, parallel = FALSE, nCpus = 4)
Error in chNet(TCGA.BRCA$X, TCGA.BRCA$group, subsampling = TRUE, R = 20,  : 
  object 'sub_iter' not found
In addition: There were 20 warnings (use warnings() to see them)
```

This built-in code without **subsampling** works great and without error.
```R
> library(chNet)
> data("TCGA.BRCA")
> result = chNet(TCGA.BRCA$X,TCGA.BRCA$group,subsampling = FALSE, R= 20,
+                lambar = 2.825, parallel = FALSE, nCpus = 4)
```
<br>

## Additional Material 
- [chNet's GitHub](https://github.com/Zhangxf-ccnu/chNet)
- [chNet's Vignette](https://github.com/Zhangxf-ccnu/chNet/blob/master/chNet_2.0.0.pdf)
- [chNet's Paper](https://academic.oup.com/bioinformatics/article/37/23/4414/6318772)
- [Supplementary Data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/37/23/10.1093_bioinformatics_btab502/2/btab502_supplementary_data.pdf?Expires=1718562724&Signature=JICple-9y7ryV5XTQVIvsoByQOYt1AIqYfz1m4qYuP6pg9DOBWWhf8OiZF-q2pqw3djPJiZ7Ep5Z1hbxeI4KRslopjKZHJyDHcGgrGP0LNesQ8R8o8FE2~pyVf6zqi0xtjPOg3qCNMNVK0-f7m4EIiVEcfcKHh2N8fmi~Ahaw2enHpe3V6nj3sV19oCsOtfmNhwRr4RqrrHgxmGeiappv17EZDbtqI2Jo1MhpZxP9P89uAAmpB81PMppYJGPm03cg7d2Ih~i~hvjYB2OQYhNJQfgHF4Yk6iuo1CAN0gsC~1nsbX4bCOvtU8wdV57gGflfyLDe8GUMxJfUJgSlyHosg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

<br>

## Reference
<a id="1">[1]</a>
Tu, Jia-Juan. (2024).
*chNet: Differential network analysis by simultaneously considering changes in gene interactions and gene expression.*
R package version 2.0.0.

<a id="2">[2]</a>
Csárdi, Gábor; Nepusz, Tamás. (2024).
*igraph: Network Analysis and Visualization in R.*
DOI: 10.5281/zenodo.7682609, R package version 2.0.3.
[igraph](https://CRAN.R-project.org/package=igraph)

<a id="3">[3]</a>
Wickham, Hadley. (2016).
*ggplot2: Elegant Graphics for Data Analysis.*
Springer-Verlag New York.
ISBN: 978-3-319-24277-4.
[ggplot2](https://ggplot2.tidyverse.org)

<a id="4">[4]</a>
Wei, Taiyun; Simko, Viliam. (2021).
*R package 'corrplot': Visualization of a Correlation Matrix.*
Version 0.92.
[corrplot](https://github.com/taiyun/corrplot)

<a id="5">[5]</a>
Gu, Zuguang. (2016, 2022).
Complex heatmap publications in *Bioinformatics* and *iMeta.*

<a id="6">[6]</a>
Gu, Zuguang; Gu, Lei; Eils, Roland; Schlesner, Matthias; Brors, Benedikt. (2014).
*circlize implements and enhances circular visualization in R.*
Bioinformatics, 30(19), pages 2811-2812.

<a id="7">[7]</a>
Kolde, Raivo. (2019).
*pheatmap: Pretty Heatmaps.*
R package version 1.0.12.
[pheatmap](https://CRAN.R-project.org/package=pheatmap)

<a id="8">[8]</a>
Neuwirth, Erich. (2022).
*RColorBrewer: ColorBrewer Palettes.*
R package version 1.1-3.
[RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer)

<a id="9">[9]</a>
Wickham, Hadley. (2007).
*Reshaping Data with the reshape Package.*
Journal of Statistical Software, 21(12), pages 1-20.
[reshape2](http://www.jstatsoft.org/v21/i12/)

<a id="10">[10]</a>
Auguie, Baptiste. (2017).
*gridExtra: Miscellaneous Functions for "Grid" Graphics.*
R package version 2.3.
[gridExtra](https://CRAN.R-project.org/package=gridExtra)

<a id="11">[11]</a>
Wu, T; Hu, E; Xu, S; Chen, M; Guo, P; Dai, Z; Feng, T; Zhou, L; Tang, L; Zhan, L; Fu, X; Liu, S; Bo, X; Yu, G. (2021).
*clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.*
The Innovation, 2(3), 100141.

<a id="12">[12]</a>
Sievert, Carson. (2020).
*Interactive Web-Based Data Visualization with R, plotly, and shiny.*
Chapman and Hall/CRC.
ISBN: 9781138331457.
[plotly](https://plotly-r.com)

<a id="13">[13]</a>
Carlson, Marc. (2023).
*org.Hs.eg.db: Genome wide annotation for Human.*
R package version 3.18.0.

<a id="14">[14]</a>
Huber, Wolfgang; Carey, V.J.; Gentleman, R.; et al. (2015).
*Orchestrating high-throughput genomic analysis with Bioconductor.*
Nature Methods, 12, pages 115-121.
[BiocGenerics](http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html)


