import pandas as pd
from gprofiler import GProfiler

def Enrichment_Analysis(tsv_file_path):
    data = pd.read_csv(tsv_file_path, sep='\t') 
    
    #Get all unique genes
    uniquegenes = pd.concat([data['regulator'], data['target']]).unique()
   
    # Initialize g:Profiler
    gp = GProfiler(return_dataframe=True)
    
    # Enrichment Analysis with Go and KEGG --- hsapiens refer to humans
    results = gp.profile(organism='hsapiens', query=uniquegenes.tolist(), 
                         sources=['GO:BP', 'GO:CC', 'GO:MF', 'KEGG'])
    
    # Filter results for significance
    p_value_results = results['p_value'] 
    significant_results = results[p_value_results < 0.05]

    # Separate GO and KEGG results
    go_results = significant_results[significant_results['source'].str.contains('GO')]
    kegg_results = significant_results[significant_results['source'] == 'KEGG']
    
    return go_results, kegg_results


# Usage
if __name__ == "__main__":
    file_path = 'network.tsv'  # Path to your network.tsv file
    go_enrichment_results, kegg_enrichment_results = Enrichment_Analysis(file_path)
    
    # Print and save the results
    print("GO Enrichment Results:")
    print(go_enrichment_results)
    go_enrichment_results.to_csv('go_enrichment_results.tsv', sep='\t', index=False)
    
    print("\nKEGG Pathway Enrichment Results:")
    print(kegg_enrichment_results)
    kegg_enrichment_results.to_csv('kegg_enrichment_results.tsv', sep='\t', index=False)