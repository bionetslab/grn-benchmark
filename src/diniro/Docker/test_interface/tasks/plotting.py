import pandas as pd
import os
from random import randrange
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import math


class Plotting:

    def plot_network(path, number_of_modules,pval_cutoff,size_threshold_min,size_threshold_max):
        print('START task plotting.py')

        ranking = pd.read_csv(os.path.join(path, 'ranked_modules.csv'))
        modules = pd.read_csv(os.path.join(path, 'modules_all.csv'))

        # plot ranking all
        plt.figure(figsize=(25, 10))
        sns.lineplot(x="TFs", y="DC-normalized", data=ranking, sort=False, marker="v")
        sns.lineplot(x="TFs", y="copula-normalized", data=ranking, sort=False, marker="o")
        sns.lineplot(x="TFs", y="Rank", data=ranking, sort=False, marker="s")

        plt.legend(title='Ranking of all Models (Named by the TFs)', loc='upper right',
                   labels=['Rank-DiffCo', 'Rank-Copula', 'Rank'])

        sns.despine(offset=0.05, trim=True);
        plt.xticks(rotation=90)
        plt.savefig(os.path.join(path, 'ranking_plot.png'), format='png')

        # plot ranking all
        plt.figure(figsize=(25, 10))
        sns.lineplot(x="TFs", y="Rank-DiffCo", data=ranking, sort=False, marker="v")
        plt.legend(title='Ranking of all Models (Named by the TFs)', loc='upper right',
                   labels=['Rank-DiffCo'])
        sns.despine(offset=0.05, trim=True);
        plt.xticks(rotation=90)
        plt.savefig(os.path.join(path, 'ranking_plot_diffco.png'), format='png')

        # plot ranking all
        plt.figure(figsize=(25, 10))
        sns.lineplot(x="TFs", y="Rank-copula", data=ranking, sort=False, marker="o")
        plt.legend(title='Ranking of all Models (Named by the TFs)', loc='upper right',
                   labels=['Rank-Copula'])

        sns.despine(offset=0.05, trim=True);
        plt.xticks(rotation=90)
        plt.savefig(os.path.join(path, 'ranking_plot_copula.png'), format='png')

        # TODO parameters
        # TODO changed see (done)
        # TODO if ties what to pick ?

        # task filter number of modules by rank
        #print(ranking)
        modules_ranked_filtered = ranking.nsmallest(int(number_of_modules), 'Rank')
        #print(modules_ranked_filtered)
        modules_ranked_filtered = modules_ranked_filtered.drop(
            ['DC-score', 'Rank-DiffCo', 'Copula-score', 'Rank-copula', 'DC-normalized', 'copula-normalized','inter_Num'], axis=1)
        #print(modules_ranked_filtered)

        modules_ranked_filtered.to_csv(os.path.join(path, 'modules_ranked_filtered.csv'), index=False)

        to_include_names = modules_ranked_filtered.values.tolist()
        print(to_include_names)

        models_list = []
        for (columnName, columnData) in modules.iteritems():
            for x in to_include_names:
                if columnData[0] == x[0]:
                    bool_series = pd.notnull(modules[columnName])
                    l = list(modules[columnName][bool_series])
                    models_list.append(l)
        #print("*************************** here ******************************")
        #print('number of modules')
        #print(len(models_list))
        #print(models_list)
        # filter by p-value :
        # TODO let the user select it
        #p_value_cutoff = 0.05
        p_vales_path = os.path.join(path, 'p_values.txt')
        with open(p_vales_path, "rb") as fp:
            pvalues = pickle.load(fp)
        #print('pvalues')
        #print(pvalues[:5])
        print('cutoff')
        print(pval_cutoff)
        pval_cutoff = float(pval_cutoff)
        #print(type(pval_cutoff))
        new_models_list = []
        for x in models_list:
            var = [[x[0], 0]]
            for y in x[1:]:
                for w in pvalues:
                    if w[0] == x[0] and w[1] == y and w[2] <= float(pval_cutoff):
                        var.append([y, str(format(w[2], "1.1e"))])
            new_models_list.append(var)
            var = []

        #print(len(new_models_list))
        #print(new_models_list)

        # TODO again send all in 1 list with tf as unique group names (done)
        interactions_1_path = os.path.join(path, 'interactions_1_p.txt')
        interactions_2_path = os.path.join(path, 'interactions_2_p.txt')

        with open(interactions_1_path, "rb") as fp:
            interactions_1 = pickle.load(fp)

        with open(interactions_2_path, "rb") as fp:
            interactions_2 = pickle.load(fp)

        z = [i for i in interactions_1 if i not in interactions_2]
        w = [i for i in interactions_2 if i not in interactions_1]
        interactions_1 = z
        interactions_2 = w

        new_models_list1 = []
        for x in new_models_list:
            lis = []
            for y in x:
                if float(y[1])<=pval_cutoff:
                    lis.append(y)
            new_models_list1.append(lis)
        print('after p-value cuttof')
        print(len(new_models_list1))
        print(new_models_list1)
        print('\n\n')



        alln = []
        alle = []
        n_combined = []
        e_combined = []
        for l in new_models_list1:
            n = {}
            nlist = []
            e = {}
            elist = []
            factor = randrange(1000)
            for i in range(len(l)):
                if i == 0:
                    n['id'] = (i + 1) * factor
                    n['label'] = l[i][0]
                    n['shape'] = 'triangleDown'
                    n['color'] = '#0000FF'
                    n['group'] = l[0][0]
                    nlist.append(n)
                    n_combined.append(n)
                    n = {}
                    # Edges
                    #e['from'] = 1 * factor
                    #e['to'] = (i + 1) * factor
                    #e['width'] = 3
                    #e['label'] = l[i][1]

                    #elist.append(e)
                    #e_combined.append(e)
                    #e = {}


                else:
                    edge_color = 'gray'
                    for x in interactions_1:
                        if x.split('_')[0] == l[0][0] and x.split('_')[1] == l[i][0]:
                            edge_color = '#008000'

                    for y in interactions_2:
                        if y.split('_')[0] == l[0][0] and y.split('_')[1] == l[i][0]:
                            edge_color = 'red'

                    if l[i][0] == l[0][0]:
                        e['from'] = 1 * factor
                        e['to'] = 1 * factor
                        e['width'] = 3
                        e['label'] = l[i][1]
                        e['color'] = edge_color
                        e['label'] = l[i][1]
                        e['arrows'] = 'to'

                        elist.append(e)
                        e_combined.append(e)
                        e = {}

                    else:
                        n['id'] = (i + 1) * factor
                        # n['label'] = '<b>' + l[i] + '</b>'
                        n['label'] = l[i][0]
                        n['shape'] = 'dot'
                        n['color'] = '#ffa500'
                        n['size'] = 20
                        n['group'] = l[0][0]
                        n['borderWidth'] = 0
                        n['borderWidthSelected'] = 2
                        nlist.append(n)
                        n_combined.append(n)
                        n = {}
                        # Edges
                        e['from'] = 1 * factor
                        e['to'] = (i + 1) * factor
                        e['arrows'] = 'to'
                        e['color'] = edge_color
                        e['length'] = 150
                        e['width'] = 3
                        e['label'] = l[i][1]

                        elist.append(e)
                        e_combined.append(e)
                        e = {}

            alln.append(nlist)
            alle.append(elist)
        # TODO this is new
        # this for plotting as one network
        intercom = []
        nod = []
        nids = []
        tfs = []
        gene_p_value_cutoff = []
        edge_p_value_cutoff = []

        #print('fixing here')
        for lis in pvalues:
            if lis[2] <= pval_cutoff:
                edge_p_value_cutoff.append(lis[0]+str('>')+lis[1])
                gene_p_value_cutoff.append(lis[0])
                gene_p_value_cutoff.append(lis[1])
                tfs.append(lis[0])
                intercom.append([lis[0], lis[1]])
        tfs = list(set(tfs))
        #print(intercom)
        #filter by size:
        count = []
        for x in intercom:
            count.append(x[0])
        print('modules sizes')
        module_pass = []
        for y in tfs:
            print(y)
            print(count.count(y))
            if count.count(y) >= size_threshold_min and count.count(y) <= size_threshold_max:
                module_pass.append(y)
        print('modules that pass the size filter')
        print(module_pass)
        print('This TFs')
        print(tfs)
        # for f in tfs:
        #     intercom.append([f, f.split('(TF)')[0]])
        #     intercom.append([f, f])

        intercom_filtred = []
        for x in intercom:
            if x[0] in module_pass:
                intercom_filtred.append(x)
        intercom = intercom_filtred

        print('intercom')
        print('\n')
        print(intercom)
        print('intercom')
        print('\n')
        gene_p_value_cutoff1 = list(set(gene_p_value_cutoff))
        edge_p_value_cutoff1 = list(set(edge_p_value_cutoff))
        print('the saved genes after p-valus')
        print(gene_p_value_cutoff1)

        gene_p_value_cutoff_path = os.path.join(path,'gene_p_value_cutoff.txt')
        with open(gene_p_value_cutoff_path, "wb") as fp:
            pickle.dump(gene_p_value_cutoff1, fp)
        #for edges
        edge_p_value_cutoff_path = os.path.join(path,'edges_p_value_cutoff.txt')
        with open(edge_p_value_cutoff_path, "wb") as fp:
            pickle.dump(edge_p_value_cutoff1, fp)

        #keep
        interactom_path = os.path.join(path,'interactom.txt')
        with open(interactom_path, "wb") as fp:
            pickle.dump(intercom, fp)

        print('************** HERE ********************')
        #print(intercom)
        for x in intercom:
            nod.append(x[0])
            nod.append(x[1])
        nodes = list(set(nod))
        for i in range(len(nodes)):
            nids.append([nodes[i],i])

        Nlist = []
        Elist = []
        nn = {}
        ee = {}
        for ww in nids:
            nn['id'] = ww[1]
            nn['label'] = ww[0]
            if ww[0] in tfs:
                nn['shape'] = 'triangleDown'
                nn['color'] = '#0000FF'
            else:
                nn['shape'] = 'circle'
                nn['color'] = '#ffdab9'
            Nlist.append(nn)
            nn = {}

        for zz in intercom:
            for y in nids:
                if zz[0] == y[0]:
                    ee['from'] = y[1]
                    f = zz[0]

                if zz[1] == y[0]:
                    ee['to'] = y[1]
                    t = zz[1]

            e_color = 'gray'
            for x in interactions_1:
                if x.split('_')[0] == f and x.split('_')[1] == t:
                    e_color = '#008000'
            for y in interactions_2:
                if y.split('_')[0] == f and y.split('_')[1] == t:
                    e_color = 'red'

            ee['color'] = e_color
            ee['width'] = 3
            ee['arrows'] = 'to'
            Elist.append(ee)
            ee = {}







        #print(Elist)
        # end of this task
        print('DONE task plotting.py \n\n\n')

        return Nlist, Elist, alln, alle, n_combined, e_combined
