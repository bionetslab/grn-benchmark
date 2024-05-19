### The parameters and settings for the simulations include:

**Number of Genes in a Module (N)**: The simulations considered modules of different sizes, specifically modules with 10, 50, or 100 genes, depending on the scenario being tested​​.  
**Proportion of Changed Correlations (γ)**: The simulations involved changing a specific proportion of the correlations within the modules to explore different degrees of network alteration. This proportion, denoted as γ, was set at different levels (0.1, 0.4, or 0.7) to represent small, medium, and large effects. In different scenarios, the changes included dropping correlations to zero or altering them by increasing or decreasing their values​​.  
**The compound symmetric correlation parameter (𝜌)**: quantifies the strength of the correlation between any two genes within the module.  
**The p-value of the p-norm difference test (PND)**: calculated using a non-parametric permutation approach, which assesses the significance of the observed differences in network structures between groups by comparing them against a distribution of differences generated under the null hypothesis.

### The parameters and settings for the testing serve as foundational structures for analyzing and visualizing the interactions between genes. such as:

**Correlation Matrix**: Measures the similarity between gene expression patterns, with each element representing the correlation coefficient between gene pairs.   
**Adjacency Matrix**: Converts correlation coefficients into a binary or weighted network structure, indicating the presence and strength of connections between genes.
**Topological Overlap Matrix (TOM)**: Enhances the adjacency matrix by considering not only direct connections between genes but also their shared connections, providing a deeper insight into the network's interconnectedness.

### Default (Suggested by Authors):

 Discomod authors suggested he paper uses ρ values of 0.3 or 0.7 in different simulation scenarios. These values represent the strength of correlation within the gene modules in the network structure.  While various values from 4 to 20 are tested, the recommended default value for general use is PND6, though the paper advises exploring other values based on specific dataset characteristics .