import gseapy
import pandas as pd
import os
import pickle
import loompy as lp
import plotly.express as px
from time import process_time
import concurrent.futures
import time
import json
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def go_function(l):
    all_genes_names = l[0]
    rank = l[1]
    output = l[2]
    library = ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018']

    try:
        # Call the Enrichr API
        enr = gseapy.enrichr(gene_list=all_genes_names, no_plot=True, organism='human', description='pathway',
                             gene_sets=library, cutoff=0.05)
        enr.run()

        # Extract results
        all_terms = enr.results

        # Filter and sort results
        filtered_terms = all_terms['Adjusted P-value'] <= 0.1
        GO = all_terms[filtered_terms]
        go = GO.sort_values('Adjusted P-value', ascending=True).head(10)

        thisisalist = list(go['Term'])
        gotermspvalues = list(go['Adjusted P-value'])
        gotermsandpvalue = []

        for i in range(len(thisisalist)):
            gotermsandpvalue.append(thisisalist[i] + "(P=" + str(format(gotermspvalues[i], "1.1e")) + ")")

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print(len(thisisalist))
        return [rank, all_genes_names[0], len(list(set(all_genes_names))), gotermsandpvalue, all_genes_names, gotermspvalues]

    except json.decoder.JSONDecodeError as e:
        logging.error("JSONDecodeError: %s", e)
        return [rank, all_genes_names[0], len(list(set(all_genes_names))), [], all_genes_names, []]
    except Exception as e:
        logging.error("An error occurred: %s", e)
        return [rank, all_genes_names[0], len(list(set(all_genes_names))), [], all_genes_names, []]

# Example usage
# Assuming l is provided as described
l = [
    ['Gene1', 'Gene2', 'Gene3'],  # Example gene list
    1,                            # Example rank
    'output'                      # Example output
]

result = go_function(l)
print(result)

class GoTerms:

    def go_terms(path,pval_cutoff,size_threshold_min,size_threshold_max):
        pval_cutoff = float(pval_cutoff)
        print('START task go_table.py')
        t1_start = process_time()

        output = os.path.join(path, 'Go_test')
        go_file = os.path.join(path, 'GO_terms.txt')

        # TODO look into the cutoff
        # TODO organism Human or Mouse

        # def go_function(all_genes_names):
        #     enr = gseapy.enrichr(gene_list=all_genes_names, no_plot=True, organism='Human', description='pathway',
        #                          gene_sets=library, cutoff=0.5,
        #                          outdir=output)
        #     all_terms = enr.results
        #     filtered_terms = all_terms['Adjusted P-value'] <= 0.5
        #     GO = all_terms[filtered_terms]
        #     return GO

        # TODO use the new one filtered.csv
        # TODO the plot from here
        modules = pd.read_csv(os.path.join(path, 'modules_all.csv'))

        filtered = pd.read_csv(os.path.join(path, 'modules_ranked_filtered.csv'))
        to_include_names = filtered.values.tolist()
        print(to_include_names)

        # plot diff expression

        def expression_plot(sample_1, module, sample):

            sample1 = sample_1[module]
            Total = sample1.sum(axis=0)
            Total = Total.to_frame(name='Expression')
            Total.div(sample_1.shape[0])
            names = [sample for x in range(Total.shape[0])]
            Total['Group'] = names
            Total.reset_index(drop=True, inplace=True)
            return Total

        sample_1_loom = os.path.join(path, 'sample_1_orgial-expression_matrix.loom')
        sample_2_loom = os.path.join(path, 'sample_2_orgial-expression_matrix.loom')

        lf = lp.connect(sample_1_loom, mode='r', validate=False)
        exprMat_sample_1 = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T

        lf = lp.connect(sample_2_loom, mode='r', validate=False)
        exprMat_sample_2 = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T

        for (columnName, columnData) in modules.iteritems():
            bool_series = pd.notnull(modules[columnName])
            a = expression_plot(exprMat_sample_1, list(modules[columnName][bool_series]), 'sampleOne')
            b = expression_plot(exprMat_sample_2, list(modules[columnName][bool_series]), 'sampleTwo')
            frames = [a, b]
            result = pd.concat(frames)

            # plotting
            title = 'Module of TF : ' + str(list(modules[columnName][bool_series])[0])
            fig = px.violin(result, y="Expression", x="Group", color="Group", box=True, points="all",
                            hover_data=result.columns, template="simple_white")
            fig.update_layout(title_text=title, title_x=0.5)
            fig.write_image(os.path.join(path, "DiffEx_" + list(modules[columnName][bool_series])[0] + ".pdf"))

        # plot

        # all_go = []
        # for x in to_include_names:
        #     for (columnName, columnData) in modules.iteritems():
        #         if columnData[0] == x[0]:
        #             bool_series = pd.notnull(modules[columnName])
        #             go = go_function(list(modules[columnName][bool_series]))
        #             go = go.sort_values('Adjusted P-value', ascending=True).head(5)
        #             thisisalist = list(go['Term'])
        #             all_go.append([int(x[1]), columnData[0], len(list(modules[columnName][bool_series])), thisisalist])

        datalist = []
        all_go = []

        # TODO let the user select it
        p_value_cutoff = 0.05
        p_vales_path = os.path.join(path, 'p_values.txt')
        with open(p_vales_path, "rb") as fp:
            pvalues = pickle.load(fp)

        for x in to_include_names:
            for (columnName, columnData) in modules.iteritems():
                if columnData[0] == x[0]:
                    bool_series = pd.notnull(modules[columnName])

                    # filter by p-value :
                    gene_list = list(modules[columnName][bool_series])
                    #print(gene_list)
                    new_gene_list = [x[0]]
                    for z in gene_list[1:]:
                        for w in pvalues:
                            if w[0] == x[0] and w[1] == z and w[2] <= pval_cutoff:
                                new_gene_list.append(z)
                    #print(new_gene_list)
                    if size_threshold_min <= len(new_gene_list)-1 <= size_threshold_max:
                        datalist.append([new_gene_list, int(x[1]), output])

        print(datalist)
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(go_function, datalist)
            for r in results:
                all_go.append(r)
        end = time.perf_counter()
        print(f'Finished in {round(end - start, 2)} s')

        with open(go_file, "wb") as fp:
            pickle.dump(all_go, fp)

        # TODO send rank also
        t1_stop = process_time()
        print('time GO terms:', t1_stop - t1_start)
        # print(all_go)
        print('DONE task go_table.py \n\n\n\n\n\n\n\n\n')
        return 1
