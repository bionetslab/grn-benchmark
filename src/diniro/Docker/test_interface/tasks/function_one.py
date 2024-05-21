import numpy as np
import pandas as pd
import networkx as nx
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import os
import pickle

class FunctionOne:

    def functionOne(path,size_thresholds):
        print('START task functionOne.py')

        # TODO get reed of this
        # Test data

        test_sample1 = [['a', 'j', 'a', 'q'], ['a', 'b', 'c', 'd'], [156, 146, 4651, 1646], [1, 1, 1, 1], [1, 1, 1, 1],
                        [164, 4641, 146, 14646]]
        test_sample1 = np.array(test_sample1).T
        test_sample1 = pd.DataFrame(test_sample1, columns=['tf', 'gene', 'weghit', 'Time', 'perc', 'score'])

        test_sample2 = [['a', 'r', 'u', 'a'], ['a', 'b', 'c', 'd'], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1],
                        [1, 1, 1, 1]]
        test_sample2 = np.array(test_sample2).T
        test_sample2 = pd.DataFrame(test_sample2, columns=['tf', 'gene', 'weghit', 'Time', 'perc', 'score'])

        # use test in case ...
        sample_1 = pd.read_csv(os.path.join(path, 'sample_grn_s1.csv'))
        sample_2 = pd.read_csv(os.path.join(path, 'sample_grn_s2.csv'))

        # Get regulons

        def regulons_extractor(sample):
            list_sample_1 = sample.values.tolist()
            tf = list(set([x[0] for x in list_sample_1]))
            regulons_s = []
            for x in tf:
                m = []
                for y in list_sample_1:
                    if y[0] == x:
                        m.append(y[1])
                regulons_s.append([x, m])
            return regulons_s

        regulons_s1 = regulons_extractor(sample_1)
        regulons_s2 = regulons_extractor(sample_2)

        common_modules = []
        for x in regulons_s1:
            for y in regulons_s2:
                if x[0] == y[0]:
                    s = x[1] + y[1]
                    common_modules.append([x[0], list(set(s))])
        print('common TFs:', [x[0] for x in common_modules])
        non_common_modules = []
        for x in regulons_s1:
            skip = 0
            for y in common_modules:
                if x[0] == y[0]:
                    skip = 1
            if skip != 1:
                non_common_modules.append(x)

        for x in regulons_s2:
            skip = 0
            for y in common_modules:
                if x[0] == y[0]:
                    skip = 1
            if skip != 1:
                non_common_modules.append(x)

        print('NON common TFs:', [x[0] for x in non_common_modules])
        all_reg = non_common_modules + common_modules
        print('ALL TFs:', [x[0] for x in all_reg])

        sample_1_filtered = sample_1.drop(columns=['Occurrence', 'occurrence_percentage'])
        sample_2_filtered = sample_2.drop(columns=['Occurrence', 'occurrence_percentage'])

        sample_1_filtered = sample_1_filtered.values.tolist()
        sample_2_filtered = sample_2_filtered.values.tolist()

        # interaction source
        sample_1_f = sample_1.drop(columns=['Weights','Occurrence', 'occurrence_percentage'])
        sample_2_f = sample_2.drop(columns=['Weights','Occurrence', 'occurrence_percentage'])

        sample_1_ff = sample_1_f.values.tolist()
        sample_2_ff = sample_2_f.values.tolist()

        only_s1 = [item for item in sample_1_ff if item not in sample_2_ff]
        only_s2 = [item for item in sample_2_ff if item not in sample_1_ff]
        both_samples = [item for item in sample_1_ff if item in sample_2_ff]

        print("interaction sources")
        print("only in sample 1", len(only_s1))
        #print(only_s1)
        print("only in sample 2",len(only_s2))
        #print(only_s2)
        print("in both ",len(both_samples))
        #print(both_samples)

        print('use this')

        only_s1_genes = []
        only_s2_genes = []
        both_samples_genes = []

        only_s1_edges = []
        only_s2_edges = []
        both_samples_edges = []

        for x in sample_1_ff:
            only_s1_edges.append(x[0]+str('>')+x[1])
            only_s1_genes.append(x[0])
            only_s1_genes.append(x[1])

        for y in sample_2_ff:
            only_s2_edges.append(y[0] + str('>')+y[1])
            only_s2_genes.append(y[0])
            only_s2_genes.append(y[1])

        only_s1_edges_f = list(set(only_s1_edges))
        only_s2_edges_f = list(set(only_s2_edges))

        only_s1_genes_f = list(set(only_s1_genes))
        only_s2_genes_f = list(set(only_s2_genes))

        # only_s1_use_edges = [item for item in only_s1_edges_f if item not in only_s2_edges_f]
        # only_s2_use_edges = [item for item in only_s2_edges_f if item not in only_s1_edges_f]
        only_s1_use_edges = only_s1_edges_f
        only_s2_use_edges = only_s2_edges_f
        both_samples_use_edges = [item for item in only_s1_edges_f if item in only_s2_edges_f]

        only_s1_use = [item for item in only_s1_genes_f if item not in only_s2_genes_f]
        only_s2_use = [item for item in only_s2_genes_f if item not in only_s1_genes_f]
        both_samples_use = [item for item in only_s1_genes_f if item in only_s2_genes_f]

        print('\n')
        print('Sample 1')
        print(only_s1_use)
        print(only_s1_use_edges)
        print('\n')
        print('Sample 2')
        print(only_s2_use)
        print(only_s2_use_edges)
        print('\n')
        print('Sample both')
        print(both_samples_use)
        print(both_samples_use_edges)
        print('\n')

        genes_sample_1_path = os.path.join(path, 'genes_sample_1.txt')
        genes_sample_2_path = os.path.join(path, 'genes_sample_2.txt')
        genes_both_path = os.path.join(path, 'genes_both.txt')

        edges_sample_1_path = os.path.join(path, 'edges_sample_1.txt')
        edges_sample_2_path = os.path.join(path, 'edges_sample_2.txt')
        edges_both_path = os.path.join(path, 'edges_both.txt')

        with open(genes_sample_1_path, "wb") as fp:
            pickle.dump(only_s1_use, fp)
        with open(genes_sample_2_path, "wb") as fp:
            pickle.dump(only_s2_use, fp)
        with open(genes_both_path, "wb") as fp:
            pickle.dump(both_samples_use, fp)

        with open(edges_sample_1_path, "wb") as fp:
            pickle.dump(only_s1_use_edges, fp)
        with open(edges_sample_2_path, "wb") as fp:
            pickle.dump(only_s2_use_edges, fp)
        with open(edges_both_path, "wb") as fp:
            pickle.dump(both_samples_use_edges, fp)

        # Done

        def network(data):
            """take data and return a network"""

            # TODO see direct graph
            G = nx.DiGraph()
            n = []
            for x in data:
                n.append(x[0])
                n.append(x[1])
            n = set(n)
            for y in n:
                G.add_node(y)
            for z in data:
                G.add_edge(z[0], z[1], weight=z[2])

            return G

        sample_1_net = network(sample_1_filtered)
        sample_2_net = network(sample_2_filtered)

        # venn diagram
        def interactionMaker(Ed):
            out = []
            for x in Ed:
                out.append(str(x[0]) + str('_') + str(x[1]))
                out.append(str(x[0]) + str('_') + str(x[0]))
            return out

        # TODO plot this in summary
        # TODO plot this after filtring
        C1_e = set(interactionMaker(sample_1_net.edges()))
        C2_e = set(interactionMaker(sample_2_net.edges()))

        interactions_1 = list(C1_e)
        interactions_2 = list(C2_e)
        interactions_1_path = os.path.join(path, 'interactions_1.txt')
        interactions_2_path = os.path.join(path, 'interactions_2.txt')
        with open(interactions_1_path, "wb") as fp:
            pickle.dump(interactions_1, fp)
        with open(interactions_2_path, "wb") as fp:
            pickle.dump(interactions_2, fp)

        venn2(subsets=[C1_e, C2_e],
              set_labels=('sample-One(' + str(len(C1_e)) + ')', 'sample-Two(' + str(len(C2_e)) + ')'),
              set_colors=('r', 'g'), alpha=0.5)
        venn2_circles(subsets=[C1_e, C2_e], linestyle='dashed', linewidth=0.6, color='k', alpha=0.5)
        plt.title('Interactions Overlaps')
        plt.savefig(os.path.join(path, 'common_interactions.png'), format='png')

        # make the 2 networks are same nodes

        missing_nodes_sample1 = (n for n in sample_2_net if n not in sample_1_net)
        sample_1_net.add_nodes_from(missing_nodes_sample1)

        missing_nodes_sample2 = (n for n in sample_1_net if n not in sample_2_net)
        sample_2_net.add_nodes_from(missing_nodes_sample2)

        nodes_order = list(sample_1_net.nodes())

        df1 = nx.to_pandas_adjacency(sample_1_net, nodelist=nodes_order, weight='weight', nonedge=0.0, dtype=float)
        df2 = nx.to_pandas_adjacency(sample_2_net, nodelist=nodes_order, weight='weight', nonedge=0.0, dtype=float)
        differential_GRN = df1 - df2
        differential_GRN = differential_GRN.abs()
        #print('differential_GRN')
        #print(differential_GRN.describe())
        # Matrix to network

        G = nx.from_pandas_adjacency(differential_GRN, create_using=nx.DiGraph)
        G.name = "Graph from pandas adjacency matrix"

        tf = []
        for x in all_reg:
            tf.append(x[0])

        sub = []
        for i, row in differential_GRN.iterrows():
            for y in all_reg:
                if i == y[0]:
                    l = list(row)
                    ll = []
                    for w in y[1]:
                        index = nodes_order.index(w)
                        ll.append(l[index])

                    sub.append([i, ll])

        # def sub_plot(d):
        #     n = []
        #     n = d[1].copy()
        #     n.append(d[0])
        #     e = [tuple([d[0], x]) for x in d[1]]
        #     G = nx.Graph()
        #     G.add_nodes_from(n)
        #     G.add_edges_from(e)
        #     nx.draw(G, pos=nx.spring_layout(G), with_labels=True, node_size=500, node_color='orange',
        #             edgecolors='black', font_size=5, edge_color='green', width=1.5)
        #     plt.show()
        #     return G
        #
        # for x in all_reg:
        #     sub_plot(x)

        # differential co-expression scores
        # TODO see if we need ascending
        dc_s = []

        for x in sub:
            dc_s.append([x[0], sum(x[1]) / len(x[1])])

        dc_s = pd.DataFrame(dc_s, columns=['TFs', 'DC-score'])
        dc_s.sort_values(by='DC-score', ascending=False)

        dc_s.to_csv(os.path.join(path, 'DC_scores.csv'), index=False)

        all_modules_l = []
        for x in all_reg:
            l = x[1].copy()
            l.insert(0, x[0])
            all_modules_l.append(l)

        modules_all = pd.DataFrame(all_modules_l)
        modules_all = modules_all.T
        #print(modules_all.count())

        print('TODO : filter by size after P-value (NOT Here))')
        # #print(df.count())
        # todrop = []
        # for i in range(len(list(modules_all.count()))):
        #     if list(modules_all.count())[i] < int(size_thresholds):
        #         todrop.append(i)
        #
        # print(todrop)
        # modules_all.drop(todrop, axis='columns', inplace=True)
        # print(modules_all.count())
        print('done')

        modules_all.to_csv(os.path.join(path, 'modules_all.csv'), index=False)
        print('DONE task functionOne.py \n\n\n')
