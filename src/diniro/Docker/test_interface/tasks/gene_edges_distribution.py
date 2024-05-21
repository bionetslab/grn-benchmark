import os
import pickle
import pandas as pd



class GeneEdgesDist:

    def update(path,pval_cutoff):
        print('New Task')
        # P values filter
        p_vales_path = os.path.join(path, 'p_values.txt')
        with open(p_vales_path, "rb") as fp:
            pvalues = pickle.load(fp)
        print('cutoff for the new task TODO: bring all from the view to here ')
        print(pval_cutoff)
        pval_cutoff = float(pval_cutoff)
        print('Here is everything after p-value')
        for x in pvalues:
            if x[2]<= pval_cutoff:
                print(x)
        new_data = []
        for w in pvalues:
            if w[2] <= pval_cutoff:
                new_data.append(w)

        #Chnage distribution after this

        interactions_1_path = os.path.join(path, 'interactions_1.txt')
        interactions_2_path = os.path.join(path, 'interactions_2.txt')

        with open(interactions_1_path, "rb") as fp:
            interactions_1 = pickle.load(fp)

        with open(interactions_2_path, "rb") as fp:
            interactions_2 = pickle.load(fp)

        interactions_1_p = []
        interactions_2_p = []

        for x in interactions_1:
            for w in new_data:
                if x.split('_')[0] == w[0] and x.split('_')[1] == w[1]:
                    interactions_1_p.append(x)
        for y in interactions_2:
            for w in new_data:
                if y.split('_')[0] == w[0] and y.split('_')[1] == w[1]:
                    interactions_2_p.append(y)

        z = [i for i in interactions_1_p if i not in interactions_2_p]
        w = [i for i in interactions_2_p if i not in interactions_1_p]
        interactions_1_p = z
        interactions_2_p = w

        interactions_1_path = os.path.join(path, 'interactions_1_p.txt')
        interactions_2_path = os.path.join(path, 'interactions_2_p.txt')
        with open(interactions_1_path, "wb") as fp:
            pickle.dump(interactions_1_p, fp)
        with open(interactions_2_path, "wb") as fp:
            pickle.dump(interactions_2_p, fp)

    def vann_diagrams(path):
        print('New Task 2')
        #TODO get interactom and compart it to the genes and edges distribution
        interactom_path = os.path.join(path, 'interactom.txt')
        with open(interactom_path, "rb") as fp:
            interactom = pickle.load(fp)

        #make inteactom looks like edges data
        interactomData = []
        for x in interactom:
            interactomData.append(x[0]+str('>')+x[1])


        edges_1_path =         interactom_path = os.path.join(path, 'edges_sample_1.txt')
        edges_2_path =         interactom_path = os.path.join(path, 'edges_sample_2.txt')

        with open(edges_1_path, "rb") as fp:
            edges_1 = pickle.load(fp)

        with open(edges_2_path, "rb") as fp:
            edges_2 = pickle.load(fp)

        # make changes to the edges distribution after the filters

        edges1Data = [i for i in edges_1 if i in interactomData]
        edges2Data = [i for i in edges_2 if i in interactomData]
        edgesBothData = [i for i in edges1Data if i in edges2Data]
        edges1Data = [i for i in edges1Data if i not in edges2Data]
        edges2Data  = [i for i in edges2Data if i not in edges1Data]
        edges1Data = [i for i in edges1Data if i not in edgesBothData]
        edges2Data  = [i for i in edges2Data if i not in edgesBothData]
        # now from edges to genes:
        gene1Data = []
        gene2Data = []
        geneBothData = []

        for x in edges1Data:
            gene1Data.append(x.split('>')[0])
            gene1Data.append(x.split('>')[1])

        for x in edges2Data:
            gene2Data.append(x.split('>')[0])
            gene2Data.append(x.split('>')[1])

        for x in edgesBothData:
            geneBothData.append(x.split('>')[0])
            geneBothData.append(x.split('>')[1])


        common = [i for i in gene1Data if i in gene2Data]
        geneBothDataExt = geneBothData + common
        gene1Data = [i for i in gene1Data if i not in geneBothDataExt]
        gene2Data = [i for i in gene2Data if i not in geneBothDataExt]
        geneBothData = geneBothDataExt
        # make things unique
        gene1Data = list(set(gene1Data))
        print(gene1Data)
        print("\n\n")
        gene2Data = list(set(gene2Data))
        print(gene2Data)
        print("\n\n")
        geneBothData = list(set(geneBothData))
        print(geneBothData)
        print("\n\n")
        edges1Data = list(set(edges1Data))
        print(edges1Data)
        print("\n\n")
        edges2Data = list(set(edges2Data))
        print(edges2Data)
        print("\n\n")
        edgesBothData = list(set(edgesBothData))
        print(edgesBothData)

        As = []
        At = []
        Bs = []
        Bt = []
        Cs = []
        Ct = []
        for x in edges1Data:
            As.append(x.split('>')[0])
            At.append(x.split('>')[1])
        for x in edges2Data:
            Bs.append(x.split('>')[0])
            Bt.append(x.split('>')[1])
        for x in edgesBothData:
            Cs.append(x.split('>')[0])
            Ct.append(x.split('>')[1])

        #mak a csv table for download
        import csv
        from itertools import zip_longest
        d = [gene1Data,gene2Data,geneBothData]
        export_data = zip_longest(*d, fillvalue='')
        with open('table.csv', 'w', encoding="ISO-8859-1", newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(("Sample I", "Sample II","Both Samples"))
            wr.writerows(export_data)
        myfile.close()
        df = pd.read_csv('table.csv')
        print(df)
        df.to_csv(os.path.join(path, 'table.csv'), index=False, header=True)

        return gene1Data,gene2Data,geneBothData,edges1Data,edges2Data,edgesBothData,As,At,Bs,Bt,Cs,Ct