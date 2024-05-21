import glob
import os
import re

import pandas as pd


class FilteringProcessing:

    def network_function(regulons, size_threshold_v):
        print('START task Filtering_processing.py: network_function')

        '''This function take regulons and filter them by wights and size (TODO: order is important)
        just by size now more then 2 coz we are using the subsampling for importance so no need for weights'''
        namesd = []
        network = []
        for x in regulons:
            namesd.append(x.name)
            network.append(x.gene2weight)

        print('Number of regulons :' + str(len(namesd)))
        for i in range(len(namesd)):
            namesd[i] = namesd[i].split(' ')
        for i in range(len(namesd)):
            namesd[i][1] = re.findall('\d+', namesd[i][1])
        for i in range(len(namesd)):
            namesd[i][1] = int(namesd[i][1][0])

        # filter by the module size
        namesd = [x for x in namesd if x[1] >= size_threshold_v]
        print('Number of regulons after size cute off:' + str(len(namesd)))

        for i in range(len(network)):
            network[i] = dict(network[i])

        t = []
        tt = []
        for i in range(len(network)):
            for x in network[i].items():
                t.append(x)
            tt.append(t)
            t = []
        data = []

        for i in range(len(namesd)):
            for j in range(len(tt[i])):
                data.append([namesd[i][0], (tt[i][j])[0], (tt[i][j])[1]])

        return data

    def open_regulon_files(name):
        regulons = []
        p_regul_p = name
        unpickled_df = pd.read_pickle(p_regul_p)
        regulons = [r.rename(r.name.replace('(+)', ' (' + str(len(r)) + 'g)')) for r in unpickled_df]
        return regulons

    # TODO threshold
    def merging_and_filtring(path, occurrence_threshold, size_thresholds):
        print('START task Filtering_processing.py: merging_and_filtring')

        all_reg = []
        size_threshold_v = 0
        #print("size threshold :",size_thresholds)
        occurrence_threshold = int(occurrence_threshold)
        print("occurrence threshold :",occurrence_threshold)

        for file in glob.glob(os.path.join(path, "*regulons.p")):
            reg = FilteringProcessing.open_regulon_files(file)
            reg_filtered = FilteringProcessing.network_function(reg, size_threshold_v)
            all_reg.append(reg_filtered)

        all_reg_list = []
        for x in all_reg:
            all_reg_list = all_reg_list + x
        # TODO look into this (section)
        if len(all_reg_list) > 0:

            pd.set_option('display.max_rows', None)
            df = pd.DataFrame(all_reg_list, columns=['TFs', 'Genes', 'Weights'])
            print('shape :', df.shape)
            #print(df)
            print('%%%%%%%%%%%%%%%%%%')
            print('Occurrence count!')
            df2 = df.groupby(["TFs", "Genes"]).size().reset_index(name="Occurrence")
            print("shape : "+ str(df2.shape))
            #print(df2)
            print('%%%%%%%%%%%%%%%%%%')
            print('Filtring by occurrence!')
            print('Occurrence threshold :', occurrence_threshold)
            df2['occurrence_percentage'] = (df2['Occurrence'] / len(all_reg)) * 100
            above_thres = df2[df2["occurrence_percentage"] >= occurrence_threshold]
            df3 = df.groupby(['TFs', 'Genes'], as_index=False)['Weights'].mean()
            int_df = pd.merge(df3, above_thres, how='inner', on=['TFs', 'Genes'])
            print("shape : "+str(int_df.shape))
            #print(int_df)
            print('%%%%%%%%%%%%%%%%%%')
            #TODO penalize large hubs
            #print('Checking for small and large hubs!')
            #print(int_df['TFs'].value_counts())
            print('%%%%%%%%%%%%%%%%%%')
            grouped_sorted = int_df.groupby(['TFs']).apply(lambda x: x.sort_values(['Weights'], ascending=False)).reset_index(drop=True)
            #print(grouped_sorted)
            #print('not here ??? after p values ??')
            #int_df1 = grouped_sorted.groupby(['TFs']).head(100)
            #print(int_df1)
            #print('checking for small and large hubs after filtring')
            #print(int_df1['TFs'].value_counts())
            print('%%%%%%%%%%%%%%%%%%')
            print("Final interactions!")
            print('shape :', grouped_sorted.shape)
            print(grouped_sorted)
            print('Here is modules size distribution')
            print(grouped_sorted['TFs'].value_counts())
            print('%%%%%%%%%%%%%%%%%%')
            print('DONE task Filtering_processing.py')
            return grouped_sorted
        else:
            print('DONE task Filtering_processing.py \n\n\n')
            return print('Messages: Empty results make it error')
