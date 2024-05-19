import pandas as pd
import anndata
import scanpy as sc
import scanet as sn
import argparse
import os
import pickle


def main(input_file_1, input_file_2, output_file,num_rep_cells,threshold) :
    
    os.chdir(output_file)

    def get_gcn_network(net,module) :
        network_ = sn.co.module_to_network(net, module = module, co_cutoff=threshold)
        return network_

    def plot_gcn_network(network,hub_genes,name = "gcn_network") :
        sn.pl.plot_gcn(network_=network,hub_genes_df=hub_genes, name=name, drug_interaction=False, smooth_edges=True)

    def plot_gcn_network_with_drugs(network,hub_genes,name = "gcn+drugs") :
        sn.pl.plot_gcn(network_=network,hub_genes_df=hub_genes, name=name, drug_interaction=True, algorithm="trustrank", smooth_edges=True)

    def get_tfs(adata, modules, Mod_,specie_) :

        return sn.grn.grn_inference(adata_processed=adata,modules_df=modules, module=Mod_, specie=specie_)

    def get_grn(network, tf, threshold, condition) :
        grn_df = network[network['Gene1'].isin(tf) | network['Gene2'].isin(tf)]
        for idx, row in grn_df.iterrows():
            if row['Gene2'] in tf:
                grn_df.at[idx, 'Gene1'], grn_df.at[idx, 'Gene2'] = grn_df.at[idx, 'Gene2'], grn_df.at[idx, 'Gene1']
        
        grn_df = grn_df.rename(columns={'Gene1': 'regulator', 'Gene2': 'target', 'Module' : 'condition', 'Weight' : 'weight'})    
        column_order = ['target', 'regulator', 'condition', 'weight']
        grn_df_final = grn_df[column_order]
        grn_df_final['condition'] = condition

        # Modify the dataframe to pass it to the plot_grn function
        grn_df = grn_df.rename(columns={'regulator': 'TF', 'target': 'TG', 'weight' : 'occurrence(pct)'})
        grn_df = grn_df.drop(columns=['condition'])
        grn_df['occurrence(pct)'] *= 100
        new_order = ['TF', 'TG', 'occurrence(pct)']
        grn_df = grn_df[new_order]

        return grn_df_final,grn_df

    def plot_grn(grn_df, threshold, name="GRN_net"):
        sn.pl.plot_grn(df=grn_df, occurrence_pct=threshold*100, name=name, regulon="all", layout="None")

    # Read the exhausted TSV file into a pandas DataFrame
    exhausted_df = pd.read_csv(input_file_1, sep='\t',index_col=0)
    exhausted_df = exhausted_df.transpose()
    exhausted_df['sample_state'] = "exhausted"

    # Read the Macrophages TSV file into a pandas DataFrame
    macrophages_df = pd.read_csv(input_file_2, sep='\t', index_col=0)
    macrophages_df = macrophages_df.transpose()
    macrophages_df['sample_state'] = "macrophage"

    # Combine the two datasets into one
    full_df = pd.concat([macrophages_df, exhausted_df], ignore_index=False)
    names = full_df.index
    names = pd.DataFrame(names, columns=['names'])
    full_df.reset_index(drop=True, inplace=True)
    full_df.columns.name = None
    full_df = full_df.sample(frac=1,random_state=777)
    full_df.reset_index(drop=True, inplace=True)

    # Get the current working directory
    current_directory = os.getcwd()

    print("Current working directory:", current_directory)


    # Drop the sample_state column and save it into a new df
    sample_state_column = full_df.pop('sample_state')
    sample_state_df = pd.DataFrame(sample_state_column, columns=['sample_state'])

    # Get the data as h5ad
    adata = anndata.AnnData(full_df)

    # Assign sample metadata DataFrame to AnnData object's observations
    adata.obs = sample_state_df

    # Compute nearest-neighbors graph
    sc.pp.neighbors(adata)

    # Compute UMAP embedding
    sc.tl.umap(adata)

    # Visualize the clusters
    sn.vz.visualization(adata, method="UMAP", color="sample_state")

    cell_annoatation="sample_state"

    # Get the representative cells
    adata_r = sn.pb.representative_cells(adata, num_rep_cells=num_rep_cells, cell_anno=cell_annoatation)

    # Get the SFTpower power that makes the network satisfy the scale-free topology
    SFTpower = sn.co.plot_powers(adata_r, network_type="unsigned", cor_method="pearson")

    # GCN inference
    net = sn.co.co_expression(adata_r, network_type="unsigned", cor_method="pearson", power=SFTpower, module_merging_threshold=0.8)


    # Number of genes per module.
    df_modules = sn.co.plot_modules(net, figsize=(8,4))

    # Identify hub genes
    hub_genes_df = sn.co.hub_genes(adata=adata_r, net=net, figsize=(6,4))

    # Compute correlation matrix
    cor = sn.co.modules_to_annotation_cor(adata_r, net, figsize=(12,6), cor_method="pearson")

    # Identify the modules of interest.
    macrophage_df = cor[cor['annotation'] == 'macrophage']
    exhausted_df = cor[cor['annotation'] == 'exhausted']
    max_cor_macrophage = macrophage_df.loc[macrophage_df['cor'].idxmax()]
    max_cor_exhausted = exhausted_df.loc[exhausted_df['cor'].idxmax()]
    M1 = max_cor_macrophage['Modules']
    M2 = max_cor_exhausted['Modules']

    var = {'net': net, 'adata_r' : adata_r}
    # Write the inferred modules into a pickle file
    with open('net.pkl', 'wb') as file:
        pickle.dump(var, file)

    # Get the network of the first module of interest
    network_ = get_gcn_network(net, module = M1)

    plot_gcn_network(network_,hub_genes_df,"gcn_network_macrophage")
    plot_gcn_network_with_drugs(network_,hub_genes_df,"gcn+drugs_macrophage")

    # Count the TFs in each module
    specie_ = 'human'
    sn.grn.regulators_count(modules_df = df_modules, specie = specie_, figsize=(8,2))

    # TFs of the first module of interest
    tfs = get_tfs(adata,df_modules,M1,specie_)

    grn_df_final, grn_df_plot = get_grn(network_,tfs, threshold,"macrophage")

    grn_df_final_macrophage = grn_df_final.reset_index(drop=True)

    plot_grn(grn_df_plot,threshold,"GRN_net_macrophages")

    # Get the network of the second module of interest
    network_ = get_gcn_network(net, module = M2)

    plot_gcn_network(network_,hub_genes_df,"gcn_network_exhausted")
    plot_gcn_network_with_drugs(network_,hub_genes_df,"gcn+drugs_exhausted")

    # TFs of the second module of interest
    tfs = get_tfs(adata,df_modules,M2,specie_)

    grn_df_final, grn_df_plot = get_grn(network_, tfs, threshold,"exhausted")

    grn_df_final_exhausted = grn_df_final.reset_index(drop=True)

    network_df = pd.concat([grn_df_final_exhausted, grn_df_final_macrophage], ignore_index=True) 

    # save the grn_df as tsv
    network_df.to_csv('network.tsv',sep='\t')

    plot_grn(grn_df_plot,threshold,"GRN_net_exhausted")

    sn.pl.plot_grn(df=grn_df_plot, occurrence_pct=threshold*100, regulon="all", name="GRN_net_with_drugs_exhausted" ,layout="None", drug_interaction="direct")

    # Get the current working directory
    current_directory = os.getcwd()
    print("You can find the grn tsv files in the following directory : " + current_directory)
    print("Gcn and grn plots can be found in outs_scan folder in the same directory")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='scanet')
    parser.add_argument('--c', type=str, help='Path to input file exhausted')
    parser.add_argument('--d', type=str, help='Path to input file macrophages')
    parser.add_argument('--o', type=str, help='Path to output directory')
    parser.add_argument('--n', type=int, default=75, help='Number of representative cells (default: 75)')
    parser.add_argument('--t', type=float, default=0.8, help='correlation_cutoff (default: 0.8)')

    # Parse the command-line arguments
    args = parser.parse_args()
    main(args.c, args.d, args.o, args.n, args.t)