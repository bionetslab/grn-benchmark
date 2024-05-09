import argparse
import pandas as pd
import os
from gprofiler import GProfiler
import gseapy as gp

# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Downstream analysis")
    parser.add_argument("--c", required=True, type=str, help="Path to the input network TSV file")
    parser.add_argument("--o", required=True, type=str, help="Directory where output files will be saved")
    return parser.parse_args()

# Load the network and split it into two networks
def load_and_split_network(input_file):
    network_df = pd.read_csv(input_file, sep='\t')
    condition_column = 'condition'
    condition_exhausted = 'exhausted'
    condition_macrophage = 'macrophage'
    exhausted_df = network_df[network_df[condition_column] == condition_exhausted].reset_index(drop=True)
    macrophage_df = network_df[network_df[condition_column] == condition_macrophage].reset_index(drop=True)
    return exhausted_df, macrophage_df

# GO enrichment analysis function
def go_enrichment(data, condition):
    targets = data['target'].dropna().tolist()
    regulators = data['regulator'].dropna().tolist()
    gene_list = list(set(targets + regulators))  # Extract unique gene identifiers
    gp = GProfiler(return_dataframe=True)
    results = gp.profile(organism='hsapiens', query=gene_list, sources=['GO:BP', 'GO:MF', 'GO:CC'])
    results['condition'] = condition
    return results

def get_all_genes(network_df):
    target_genes = network_df['target'].unique().tolist()
    regulator_genes = network_df['regulator'].unique().tolist()
    all_genes = list(set(target_genes + regulator_genes))
    return all_genes

# KEGG enrichment analysis function
def kegg_analysis(gene_list, output_dir, condition_name):
    kegg_enrichment = gp.enrichr(
        gene_list=gene_list,
        gene_sets='KEGG_2021_Human',
        outdir=output_dir,
        no_plot=True
    )
    kegg_enrichment.results['condition'] = condition_name
    return kegg_enrichment.results

# Main function
def main():
    args = parse_args()
    os.makedirs(args.o, exist_ok=True)
    os.chdir(args.o)

    exhausted_df, macrophage_df = load_and_split_network(args.c)

    # GO Analysis
    exhausted_results = go_enrichment(exhausted_df, 'exhausted')
    macrophage_results = go_enrichment(macrophage_df, 'macrophage')
    concatenated_results = pd.concat([exhausted_results, macrophage_results])
    concatenated_results.to_csv("Go_enriched_terms.tsv", sep="\t", index=False)

    concatenated_results['ratio'] = concatenated_results['intersection_size'] / concatenated_results['term_size']
    top_go_enriched_terms = concatenated_results.groupby('condition').apply(
        lambda x: x.nlargest(5, 'ratio')).reset_index(drop=True)
    top_go_enriched_terms.to_csv("Top_go_enriched_terms.tsv", sep="\t", index=False)

    # KEGG Analysis
    genes_exhausted = get_all_genes(exhausted_df)
    genes_macrophage = get_all_genes(macrophage_df)
    results_exhausted = kegg_analysis(genes_exhausted, args.o, 'exhausted')
    results_macrophage = kegg_analysis(genes_macrophage, args.o, 'macrophage')
    combined_results = pd.concat([results_exhausted, results_macrophage], ignore_index=True)
    combined_results.to_csv("KEGG_enriched_pathways.tsv", sep="\t", index=False)

    top_pathways = combined_results.groupby('condition').apply(
        lambda x: x.nlargest(5, 'Combined Score')
    ).reset_index(drop=True)
    top_pathways.to_csv("Top_KEGG_enriched_pathways.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
