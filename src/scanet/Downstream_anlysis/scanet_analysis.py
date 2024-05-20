import scanet as sn
import pickle
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Modules analysis")
    parser.add_argument('--c', type=str, help='Path to the pickle file')
    return parser.parse_args()

def analyze_module(M, adata_r, module_membership, cor):
    sn.co.module_to_annotation_cor(cor, module=M, figsize=(8, 4))
    sn.co.Module_Genes_Avrage_Expr(module=M, adata=adata_r, figsize=(12, 6))
    sn.co.Module_Activity(module=M, plot_type="box", df=module_membership, figsize=(12, 6))
    sn.co.Module_Activity(module=M, plot_type="violin", df=module_membership, figsize=(12, 6))

def main():
    args = parse_args()
    pickle_file_path = args.c

    if not os.path.isfile(pickle_file_path):
        print(f"The file {pickle_file_path} does not exist.")
        return

    with open(pickle_file_path, 'rb') as file:
        pickle_file = pickle.load(file)

    net = pickle_file['net']
    adata_r = pickle_file['adata_r']

    if not adata_r.obs['__SCANclusters__'].dtype == 'category':
        adata_r.obs['__SCANclusters__'] = adata_r.obs['__SCANclusters__'].astype('category')

    sn.co.plot_dendrogram(net, fig_size=(800, 400))
    sn.co.plot_eigengene_network(net, fig_size=(800, 400))  # Heatmap showing pairwise correlations between module eigengenes
    modules = sn.co.plot_modules(net, figsize=(8, 4))  # Number of genes per module
    print(modules)
    hub_genes = sn.co.hub_genes(adata=adata_r, net=net, figsize=(6, 4))  # Hub genes
    print(hub_genes)
    cor = sn.co.modules_to_annotation_cor(adata_r, net, figsize=(12, 6), cor_method="pearson")  # Correlation between the modules and the two conditions
    print(cor)
    module_membership = sn.co.plot_module_membership(net, adata=adata_r, figsize=(20, 20))

    macrophage_df = cor[cor['annotation'] == 'macrophage']
    exhausted_df = cor[cor['annotation'] == 'exhausted']
    max_cor_macrophage = macrophage_df.loc[macrophage_df['cor'].idxmax()]
    max_cor_exhausted = exhausted_df.loc[exhausted_df['cor'].idxmax()]
    M_macrophage = max_cor_macrophage['Modules']
    M_exhausted = max_cor_exhausted['Modules']
    
    # Analyze the module associated with the macrophage pretreatment condition
    analyze_module(M_macrophage, adata_r, module_membership, cor)
    
    # Analyze the module associated with the exhausted condition
    analyze_module(M_exhausted, adata_r, module_membership, cor)

if __name__ == "__main__":
    main()
