INTRODUCTION TO THE TOOL
Background
Both differential expression (DE) and differential co-expression (DC) analyses are appreciated as useful tools in understanding gene regulation related to complex diseases. The performance of integrating DE and DC, however, remains unexplored.

Results
In this study, we proposed a novel analytical approach called DECODE (Differential Co-expression and Differential Expression) to integrate DC and DE analyses of gene expression data. DECODE allows one to study the combined features of DC and DE of each transcript between two conditions. By incorporating information of the dependency between DC and DE variables, two optimal thresholds for defining substantial change in expression and co-expression are systematically defined for each gene based on chi-square maximization. By using these thresholds, genes can be categorized into four groups with either high or low DC and DE characteristics. In this study, DECODE was applied to a large breast cancer microarray data set consisted of two thousand tumor samples. By identifying genes with high DE and high DC, we demonstrated that DECODE could improve the detection of some functional gene sets such as those related to immune system, metastasis, lipid and glucose metabolism. Further investigation on the identified genes and the associated functional pathways would provide an additional level of understanding of complex disease mechanism.

Conclusions
By complementing the recent DC and the traditional DE analyses, DECODE is a valuable methodology for investigating biological functions of genes exhibiting disease-associated DE and DC combined characteristics, which may not be easily revealed through DC or DE approach alone.

DECODE is available at the Comprehensive R Archive Network (CRAN): http://cran.r-project.org/web/packages/decode/index.html.

LINK TO THE PUBLICATION : https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0582-4

METHODS
DECODE consists of four steps: (1) calculating differential expression (DE), (2) calculating differential co-expression (DC), (3) selecting thresholds to define high or low values of DC and DE variables based on chi-square maximization, and statistically evaluating partitions divided by the thresholds, (4) comparing functional relevance of genes categorized into the partitions of high DC, high DE, or both. Figure 1 illustrates the overview of the analytical framework. Details are described in the following sections.

Package description
Integrated differential expression (DE) and differential co-expression (DC) analysis on gene expression data based on DECODE (DifferEntial CO-expression and Differential Expression) algorithm. Given a set of gene expression data and functional gene set data, the program will return a table summary of the selected gene sets with high differential co-expression and high differential expression (HDC-HDE).

INPUT PARAMETERS FOR THE TOOL
(1) gene expression data

(2) functional gene set data

Output: Table summary for the selected HDC-HDE gene sets, 'out_summary.txt'

Data format for gene expression data (Columns are tab-separated):

Column 1: Official gene symbol

Column 2: Probe ID

Starting from column 3: Expression for different samples

Row 1 (starting from column 3): Sample class ("1" indicates control group; "2" indicates case group)

Row 2: Sample id

Starting from row 3: Expression for different genes

Example:

geneName probeID 2 2 2 1 1 1

       -          Case1    Case2    Case3   Control1  Control2  Control3
7A5 ILMN_1762337 5.12621 5.19419 5.06645 5.40649 5.51259 5.387

A1BG ILMN_2055271 5.63504 5.68533 5.66251 5.37466 5.43955 5.50973

A1CF ILMN_2383229 5.41543 5.58543 5.43239 5.49634 5.62685 5.36962

A26C3 ILMN_1653355 5.56713 5.5547 5.59547 5.46895 5.49622 5.50094

A2BP1 ILMN_1814316 5.23016 5.33808 5.31413 5.30586 5.40108 5.31855

A2M ILMN_1745607 7.65332 6.56431 8.20163 9.19837 9.04295 10.1448

A2ML1 ILMN_2136495 5.53532 5.93801 5.33728 5.36676 5.79942 5.13974

A3GALT2 ILMN_1668111 5.18578 5.35207 5.30554 5.26107 5.26536 5.28932

Data format for functional gene set data (Columns are tab-separated):

Column 1: Functional gene set name

Column 2: Other description such as gene set id

Starting from column 3: Official gene symbols for the functional gene set

Example:

B cell activation GO\GO:0042113 AKAP17A ZAP70 PFDN1 ...

apoptotic signaling pathway GO\GO:0097190 ITPR1 PTH DNAJC10 HINT1 ...

To load the package

install.packages(decode)
install.packages(tidyverse)
install.packages(ggplot2)

Initiate Execution

Example 1:
Running a larger set of gene expression data with 1400 genes. It will take ~16 minutes to complete. (Computer used: An Intel Core i7-4600 processor, 2.69 GHz, 8 GB RAM)

path = system.file('extdata', package='decode')
geneSetInputFile = file.path(path, "geneSet.txt")
geneExpressionFile = file.path(path, "Expression_data_1400genes.txt")
runDecode(geneSetInputFile, geneExpressionFile)
Example 2:
A smaller set of gene expression data with 50 genes to satisfy CRAN’s submission requirement. (No results will be generated)

path = system.file('extdata', package='decode')
geneSetInputFile = file.path(path, "geneSet.txt")
geneExpressionFile = file.path(path, "Expression_data_50genes.txt")
runDecode(geneSetInputFile, geneExpressionFile)
[1] "Reading gene expression data..."
[1] "Calculating t-statistics..."
[1] "Calculating pairwise correlation for normal states..."
[1] "Calculating pairwise correlation for disease states..."
[1] "Calculating differential co-expression measures ..."
[1] "Reading functional gene set data"
[1] "Identifying optimal thresholds for genes"
[1] "Gene id: 1"
[1] "Gene id: 2"
[1] "Gene id: 3"
[1] "Gene id: 4"
[1] "Gene id: 5"
[1] "Gene id: 6"
[1] "Gene id: 7"
[1] "Gene id: 8"
[1] "Gene id: 9"
[1] "Gene id: 10"
[1] Gene id: 11"
[1] "Gene id: 12"
[1] "Gene id: 13"
[1] "Gene id: 14"
[1] "Gene id: 15"
[1] "Gene id: 16"
[1] "Gene id: 17"
[1] "Gene id: 18"
[1] "Gene id: 19"
[1] "Gene id: 20"
[1] "Gene id: 21"
[1] "Gene id: 22"
[1] "Gene id: 23"
[1] "Gene id: 24"
[1] "Gene id: 25"
[1] "Gene id: 26"
[1] "Gene id: 27"
[1] "Gene id: 28"
[1] "Gene id: 29"
[1] "Gene id: 30"
[1] "Gene id: 31"
[1] "Gene id: 32"
[1] "Gene id: 33"
[1] "Gene id: 34"
[1] "Gene id: 35"
[1] "Gene id: 36"
[1] "Gene id: 37"
[1] "Gene id: 38"
[1] "Gene id: 39"
[1] "Gene id: 40"
[1] "Gene id: 41"
[1] "Gene id: 42"
[1] "Gene id: 43"
[1] "Gene id: 44"
[1] "Gene id: 45"
[1] "Gene id: 46"
[1] "Gene id: 47"
[1] "Gene id: 48"
[1] "Gene id: 49"
[1] "Gene id: 50"
[1] "Identifying best associated functional gene set for each gene..."
[1] "Gene id: 1"
[1] "Gene id: 2"
[1] "Gene id: 3"
[1] "Gene id: 4"
[1] "Gene id: 5"
[1] "Gene id: 6"
[1] "Gene id: 7"
[1] "Gene id: 8"
[1] "Gene id: 9"
[1] "Gene id: 10"
[1] "Gene id: 11"
[1] "Gene id: 12"
[1] "Gene id: 13"
[1] "Gene id: 14"
[1] "Gene id: 15"
[1] "Gene id: 16"
[1] "Gene id: 17"
[1] "Gene id: 18"
[1] "Gene id: 19"
[1] "Gene id: 20"
[1] "Gene id: 21"
[1] "Gene id: 22"
[1] "Gene id: 23"
[1] "Gene id: 24"
[1] "Gene id: 25"
[1] "Gene id: 26"
[1] "Gene id: 27"
[1] "Gene id: 28"
[1] "Gene id: 29"
[1] "Gene id: 30"
[1] "Gene id: 31"
[1] "Gene id: 32"
[1] "Gene id: 33"
[1] "Gene id: 34"
[1] "Gene id: 35"
[1] "Gene id: 36"
[1] "Gene id: 37"
[1] "Gene id: 38"
[1] "Gene id: 39"
[1] "Gene id: 40"
[1] "Gene id: 41"
[1] "Gene id: 42"
[1] "Gene id: 43"
[1] "Gene id: 44"
[1] "Gene id: 45"
[1] "Gene id: 46"
[1] "Gene id: 47"
[1] "Gene id: 48"
[1] "Gene id: 49"
[1] "Gene id: 50"
[1] "Processing raw results..."
[1] "Summarizing functional gene set results..."
[1] "Done. Result is saved in out_summary.txt"`
####OUTPUT

We receive an output in the text format stating the best associated gene set with highest minimum Functional Information.

###Conclusion

We presented a novel method named DECODE as a mean to integrate the DC and DE analysis. DECODE provides an analytic framework for studying different DC and DE characteristics of the genes. By incorporating dependency between DC and DE, high or low values of the DC and DE variables are systematically defined by selecting optimal thresholds that maximize the chi-square value. In using the optimal thresholds, genes can be divided into partitions with different DC and DE characteristics. The statistical significance of a gene partition can further be evaluated by residual test. Noteworthy, since the identified gene partitions at this stage are not constrained or depended on any predefined functional modules or pathways, they provide the opportunities for the discovery of novel disease related genes.

DECODE is useful for investigating whether the functional information of an identified gene partition using the combining DC and DE criteria is higher than that using individual DC or DE criteria alone. In other words, it may generate critical novel biological insights which may not be easily obtained using individual DC or DE approach. In applying DECODE to the breast cancer data, we demonstrated that it can improve the detection of some immune system, metastasis, lipid and glucose metabolism related gene sets using high DE and high DC criteria. Further investigation on the identified gene partitions and the associated functional pathways provides a more systematic understanding of complex disease mechanism, which in turn yields useful insights in the development of new therapeutic strategies for the disease. In conclusion, in complementing the DC and DE analysis, DECODE is a valuable methodology in identifying functional gene sets exhibiting certain combination of DE and DC characteristics, which serves as a new tool for future gene expression studies.

References
[1] https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0582-4