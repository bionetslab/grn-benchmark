import numpy as np
import pandas as pd
from tqdm import tqdm

from scipy.stats import ks_2samp


class GeneExpressionAnalyzer:
    TIES_AVERAGE = "average"
    TIES_MIN = "min"
    TIES_MAX = "max"
    TIES_DENSE = "dense"
    TIES_ORIGINAL = "ordinal"

    def __init__(self, empirical_copula):
        """
        Initializes the GeneExpressionAnalyzer with an instance of EmpiricalCopula.

        Args:
            empirical_copula (EmpiricalCopula): An instance of EmpiricalCopula used to compute pseudo-observations
                                                 and empirical copulas.
        """
        self.empirical_copula = empirical_copula

    def compute_dc_copula_network(
        self,
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        ties_method: str = "average",
        smoothing: str = "none",
        ks_stat_method: str = "asymp",
    ) -> pd.DataFrame:
        """
        Computes a network of differential coexpression scores for gene pairs across two conditions.
        The score for each gene pair is calculated using the Kolmogorov-Smirnov distance between their empirical
        copulas, representing the degree of differential coexpression. This network is used to assess the similarity
        in gene expression distributions between 'tumor' and 'normal' phenotypes for example.
        """
        # Extract gene names from the first column
        gene_names = df1.iloc[:, 0].values
        assert np.array_equal(
            gene_names, df2.iloc[:, 0].values
        ), "Gene lists must match!"

        # Extracting numeric data from the dataframes, assuming the first column is the header
        data1 = df1.iloc[:, 1:].values
        data2 = df2.iloc[:, 1:].values
        n_genes = min(data1.shape[0], data2.shape[0])

        # Print dataset summary
        print(f"Starting DC Copula coexpression calculation:")
        print(f" - Number of gene pairs to be analyzed: {n_genes * (n_genes - 1) // 2}")
        print(f" - Ties method: {ties_method}")
        print(f" - Smoothing technique: {smoothing}")
        print(f" - KS statistic mode: {ks_stat_method}")

        # Create a tqdm progress bar
        pbar = tqdm(
            total=(n_genes * (n_genes - 1)) // 2,
            desc="Computing distances",
            unit="pair",
        )

        # Initialize an empty list for the network
        rows_list = []

        # Loop over all unique pairs of genes, calculated as the number of combinations of n_genes taken 2 at a time
        for k in range(n_genes * (n_genes - 1) // 2):
            # Calculate the first index 'i' of the gene pair
            # This formula ensures 'i' progresses through each gene until the second to last gene
            i = int((2 * n_genes - 1 - np.sqrt((2 * n_genes - 1) ** 2 - 8 * k)) / 2)
            # Calculate the second index 'j' of the gene pair
            # This formula ensures 'j' is always greater than 'i', thus avoiding duplicate pairs and self-pairs
            j = (
                k
                + i
                + 1
                - n_genes * (n_genes - 1) // 2
                + (n_genes - i) * ((n_genes - i) - 1) // 2
            )

            # Stack the gene expression data for genes i and j from the first dataset, transposed to match samples by rows
            gene_pair_data1 = np.vstack((data1[i, :], data1[j, :])).T
            # Repeat the stacking for the second dataset
            gene_pair_data2 = np.vstack((data2[i, :], data2[j, :])).T

            # Calculate pseudo-observations for the gene pair in the first condition
            u1 = self.empirical_copula.pseudo_observations(gene_pair_data1, ties_method)
            # Calculate pseudo-observations for the gene pair in the second condition
            u2 = self.empirical_copula.pseudo_observations(gene_pair_data2, ties_method)

            # Compute the empirical copula for the first condition
            ec1 = self.empirical_copula.empirical_copula(
                u1, gene_pair_data1, ties_method, smoothing
            )
            # Compute the empirical copula for the second condition
            ec2 = self.empirical_copula.empirical_copula(
                u2, gene_pair_data2, ties_method, smoothing
            )

            # Calculate the Kolmogorov-Smirnov statistic between the two empirical copulas
            # This statistic quantifies the difference in joint distribution of expression levels between two conditions
            ks_stat, _ = ks_2samp(ec1, ec2, method=ks_stat_method)

            # Append the computed data as a dictionary to a list; each dictionary represents a connection (edge)
            # in the network, labeled with the indices of the genes, the condition, and the weight (KS statistic)
            rows_list.append(
                {
                    "Target": gene_names[j],
                    "Regulator": gene_names[i],
                    "Condition": "Diff Co-Exp between both Condition",
                    "Weight": ks_stat,
                }
            )

            pbar.update(1)  # Update the progress bar after each iteration

        pbar.close()  # Close the progress bar when done
        # Creating a DataFrame from the list of rows
        network_df = pd.DataFrame(rows_list)
        return network_df
