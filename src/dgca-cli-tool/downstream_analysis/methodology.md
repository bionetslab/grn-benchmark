## Methodology

This section describes the methodology implemented by the CLI tool to perform Differential Gene Correlation Analysis using the DGCA package. The process is structured into several steps, starting from the initial data input to the final output generation, as detailed in Algorithm 1. 

### Algorithm 1: Workflow of the DGCA Command-Line Tool

#### Requirements:
- Two gene expression files of distinct conditions

#### Output:
- `network.tsv` file containing the differential gene correlation analysis results. A detailed mathematical explanation of the `network.tsv` file's weights can be found in the [README.md](../README.md) file of my CLI tool’s repository.

#### Steps:

1. **Input File:**
   - Takes gene expression data of two conditions as input.

2. **Prepare Inputs for DGCA Package:**
   - Construct a matrix of gene expression values.
   - Create a design matrix specifying the conditions associated with samples.
   - Specify the conditions for comparison.

3. **Execute DGCA:**
   - Run the DGCA analysis using the prepared input matrices from step 2 to generate results.

4. **Format Conversion of Results:**
   - Convert the raw results from DGCA into the required format for the project.

5. **Generate Output File:**
   - Output the final results into a `network.tsv` file which contains the network of differentially correlated genes as computed by the DGCA.

6. **Return Output:**
   - The process returns `network.tsv` containing the results of the differential gene correlation analysis.


