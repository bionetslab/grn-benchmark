import argparse
import os
import scanet as sn
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Paper reproduction")
    parser.add_argument("--d", type=str, required=True, help="Path to the h5ad data file")
    parser.add_argument('--o', type=str, help='Path to output directory')
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()
    data_path = args.d
    output_dir = args.o
    os.chdir(output_dir)
    adata = sn.pp.read_h5ad(data_path,pr_process='Yes')

    # Extract the most 4000 highly expressed genes
    adata = sn.pp.extract_highly_variable_genes(adata, n_top_genes=4000, plot=True)
    adata = adata[:, adata.var.highly_variable]

    # Convert 'diet' and 'louvain' columns to strings
    adata.obs['diet'] = adata.obs['diet'].astype(str)
    adata.obs['louvain_anno'] = adata.obs['louvain_anno'].astype(str)
    # Create a new column by combining 'diet' and 'louvain'
    adata.obs['Small intestine'] = adata.obs['diet'] + '_' + adata.obs['louvain_anno']

    # Visualization
    cell_annoatation="Small intestine"
    sn.vz.visualization(adata, method="UMAP", color=cell_annoatation)

    adata_r = sn.pb.representative_cells(adata, num_rep_cells=45, cell_anno="Small intestine")
    #SFTpower = sn.co.plot_powers(adata_r, network_type="unsigned", cor_method="pearson")
    net = sn.co.co_expression(adata_r, network_type="unsigned", cor_method="pearson", power=3,module_merging_threshold=0.8)
    df_modules = sn.co.plot_modules(net, figsize=(8,4))
    cor = sn.co.modules_to_annotation_cor(adata_r, net, figsize=(12,6), cor_method="pearson")
    # Get the module of intereset
    HFD_ISC_df = cor[cor['annotation'] == 'HFD_ISC']
    max_cor_HFD_ISC = HFD_ISC_df.loc[HFD_ISC_df['cor'].idxmax()]
    M = max_cor_HFD_ISC['Modules']
    hub_genes_df = sn.co.hub_genes(adata=adata_r, net=net, figsize=(6,4))
    network_ = sn.co.module_to_network(net, module = M, co_cutoff=0.6)
    sn.pl.plot_gcn(network_=network_,hub_genes_df=hub_genes_df, name="gcn+drugs_HFD_ISC", drug_interaction=True, algorithm="trustrank", smooth_edges=True)
    specie_ = 'mouse'
    sn.grn.regulators_count(modules_df = df_modules, specie = specie_, figsize=(8,2))
    tfs = sn.grn.grn_inference(adata_processed=adata,modules_df=df_modules, module=M, specie=specie_)
    grn_df = network_[network_['Gene1'].isin(tfs) | network_['Gene2'].isin(tfs)]
    for idx, row in grn_df.iterrows():
        if row['Gene2'] in tfs:
            grn_df.at[idx, 'Gene1'], grn_df.at[idx, 'Gene2'] = grn_df.at[idx, 'Gene2'], grn_df.at[idx, 'Gene1']
    grn_df = grn_df.rename(columns={'Gene1': 'regulator', 'Gene2': 'target', 'Module' : 'condition', 'Weight' : 'weight'})
    threshold =  0.6
    column_order = ['target', 'regulator', 'condition', 'weight']
    grn_df_final = grn_df[column_order]
    grn_df_final['condition'] = 'HFD_ISC'
    grn_df_final = grn_df_final.reset_index(drop=True)
    grn_df_final.to_csv('grn_df_HFD_ISC.tsv',sep='\t')
    grn_df = grn_df.rename(columns={'regulator': 'TF', 'target': 'TG', 'weight' : 'occurrence(pct)'})
    grn_df = grn_df.drop(columns=['condition'])
    grn_df['occurrence(pct)'] *= 100
    new_order = ['TF', 'TG', 'occurrence(pct)']
    grn_df = grn_df[new_order]
    sn.pl.plot_grn(df=grn_df, occurrence_pct=threshold*100, name="GRN_net_HFD_ISC", regulon="all", layout="None")
