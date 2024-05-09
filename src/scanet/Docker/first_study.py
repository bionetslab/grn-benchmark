import argparse
import scanet as sn
import pandas as pd
import matplotlib.pyplot as plt
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Paper reproduction")
    parser.add_argument("--d", type=str, required=True, help="Path to the h5ad data file")
    parser.add_argument("--c", type=str, required=True, help="Path to the modules DataFrame CSV file")
    parser.add_argument('--o', type=str, help='Path to output directory')
    return parser.parse_args()

def plot_no_genes_per_module(modules_df):
    sorted_modules_df = modules_df.sort_values(by='Modules')
    # Extract the sorted modules and corresponding number of genes
    sorted_modules = sorted_modules_df['Modules'].tolist()
    corresponding_genes = sorted_modules_df['Number of genes'].tolist()
    # Convert Modules column to string to ensure correct sorting
    sorted_modules_df['Modules'] = sorted_modules_df['Modules'].astype(str)   
    def custom_sort(module):
        return int(module[1:])
    
    sorted_modules_df = sorted_modules_df.iloc[sorted_modules_df['Modules'].map(custom_sort).argsort()]
    sorted_modules = sorted_modules_df['Modules'].tolist()
    corresponding_genes = sorted_modules_df['Number of genes'].tolist()
    modules = sorted_modules
    num_genes = corresponding_genes
    # Create the bar plot
    plt.figure(figsize=(15, 8)) 
    plt.bar(modules, num_genes, color='b', edgecolor='black')

    for i, num in enumerate(num_genes):
        plt.text(i, num + 5, str(num), ha='center', va='bottom')

    plt.xlabel("Modules")
    plt.ylabel("Number of genes")
    plt.title("Number of genes per Module")
    plt.grid(axis='y', linestyle=' ', alpha=0.7)
    plt.tight_layout()
    plt.xticks(rotation=45)
    plt.show()

def gcn_and_grn_inference(net,module,modules_df,threshold,adata_r,condition):
    network_ = sn.co.module_to_network(net, module =module, co_cutoff=threshold)
    hub_genes_df = sn.co.hub_genes(adata=adata_r, net=net, figsize=(6,4))
    sn.pl.plot_gcn(network_=network_,hub_genes_df=hub_genes_df, name=f"gcn_network_{condition}", drug_interaction=False, smooth_edges=True)
    tfs = sn.grn.grn_inference(adata_processed=classical_monocytes_data,modules_df=modules_df, module=module, specie="human")
    grn_df = network_[network_['Gene1'].isin(tfs) | network_['Gene2'].isin(tfs)]
    for idx, row in grn_df.iterrows():
        if row['Gene2'] in tfs:
            grn_df.at[idx, 'Gene1'], grn_df.at[idx, 'Gene2'] = grn_df.at[idx, 'Gene2'], grn_df.at[idx, 'Gene1']
    grn_df = grn_df.rename(columns={'Gene1': 'regulator', 'Gene2': 'target', 'Module' : 'condition', 'Weight' : 'weight'})    
    column_order = ['target', 'regulator', 'condition', 'weight']
    grn_df_final = grn_df[column_order]
    grn_df_final['condition']=condition
    grn_df_final = grn_df_final.reset_index(drop=True)
    grn_df_final.to_csv(f'grn_df_{condition}.tsv',sep='\t')
    grn_df = grn_df.rename(columns={'regulator': 'TF', 'target': 'TG', 'weight' : 'occurrence(pct)'})
    grn_df = grn_df.drop(columns=['condition'])
    grn_df['occurrence(pct)'] *= 100
    order = ['TF', 'TG', 'occurrence(pct)']
    grn_df = grn_df[order]
    sn.pl.plot_grn(df=grn_df, occurrence_pct=threshold*100, name=f"GRN_net_{condition}", regulon="all", layout="None")


if __name__ == "__main__":

    args = parse_arguments()
    data_path = args.d
    path_to_modules_df = args.c
    output_dir = args.o

    os.chdir(output_dir)

    adata = sn.pp.read_h5ad(data_path)
    # Get the genes names
    adata.var = adata.var.set_index('features')

    # Visualization
    cell_annoatation="cell_types"
    sn.vz.visualization(adata, method="UMAP", color=cell_annoatation)

    # Get only classical_monocytes cells
    classical_monocytes_indices = adata.obs['cell_types'] == "Classical monocytes"
    classical_monocytes_data = adata[classical_monocytes_indices, :]
    classical_monocytes_data.obs['outcome'] = classical_monocytes_data.obs['outcome'].replace({'deceased': 'Deceased'})
    classical_monocytes_data.obs['outcome'] = classical_monocytes_data.obs['outcome'].replace({'NA': 'Healthy'})
    classical_monocytes_data.obs['outcome'] = classical_monocytes_data.obs['outcome'].replace({'survived': 'Survived'})
    cell_annoatation="outcome"
    sn.vz.visualization(classical_monocytes_data, method="UMAP", color=cell_annoatation)

    # Read the provided modules by the author
    modules_df = pd.read_csv(path_to_modules_df)
    modules_df = modules_df.rename(columns={"Genes": "genes"})
    modules_df.head()

    plot_no_genes_per_module(modules_df)

    # Create my own modules_df, which will be used in gcn and grn inference. The modules_df provided by the author could't be used, as the parameter setting is missing.
    adata_r = sn.pb.representative_cells(classical_monocytes_data, num_rep_cells=150, cell_anno="outcome")
    SFTpower = sn.co.plot_powers(adata_r, network_type="unsigned", cor_method="pearson")
    net = sn.co.co_expression(adata_r, network_type="unsigned", cor_method="pearson", power=SFTpower,module_merging_threshold=0.8)
    My_modules_df = sn.co.plot_modules(net, figsize=(8,4))

    cor = sn.co.modules_to_annotation_cor(adata_r, net, figsize=(12,6), cor_method="pearson")
    # Get modules of interest
    deaceased_df = cor[cor['annotation'] == 'Deceased']
    survived_df = cor[cor['annotation'] == 'Survived']
    max_cor_deaceased = deaceased_df.loc[deaceased_df['cor'].idxmax()]
    max_cor_survived = survived_df.loc[survived_df['cor'].idxmax()]
    M1 = max_cor_deaceased['Modules']
    M2 = max_cor_survived['Modules']

    specie_ = 'human'
    sn.grn.regulators_count(modules_df = My_modules_df, specie = specie_, figsize=(8,2))
    # Gcn and grn inference of the first module of interest
    gcn_and_grn_inference(net,M1,My_modules_df,threshold=0.6,adata_r=adata_r,condition="deceased")
    # Gcn and grn inference of the second module of interest
    gcn_and_grn_inference(net,M2,My_modules_df,threshold=0.85,adata_r=adata_r,condition="survived")
