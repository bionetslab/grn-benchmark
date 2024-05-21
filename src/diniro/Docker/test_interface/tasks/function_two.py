import pickle
import loompy as lp
import pandas as pd
import os
import scipy
from time import process_time
import concurrent.futures
import time
from sklearn import preprocessing
import numpy as np
from scipy.stats import rankdata, ks_2samp

import numpy as np
import pandas as pd
from scipy.stats import rankdata, ks_2samp

def empirical_copula(data):
    ranks = np.apply_along_axis(rankdata, 0, data) / (len(data) + 1)
    return ranks

def calculate_empirical_cdf(data):
    """
    Calculate the empirical CDF for copula data.
    """
    n = len(data)
    ecdf = np.arange(1, n + 1) / n
    return ecdf

def copulas(input_list):
    tf = input_list[0][0]
    l = input_list[0][1:]
    exprMat_sample_1 = input_list[1]
    exprMat_sample_2 = input_list[2]
    var = []
    inter_dis = []
    
    for x in l:
        df1 = exprMat_sample_1[[tf, x]]
        df2 = exprMat_sample_2[[tf, x]]
        
        copula_df1 = empirical_copula(df1.values)
        copula_df2 = empirical_copula(df2.values)
        
        ec1 = calculate_empirical_cdf(copula_df1[:, 1])
        ec2 = calculate_empirical_cdf(copula_df2[:, 1])
        
        ks_statistic, _ = ks_2samp(ec1, ec2)
        Dis = 1 - ks_statistic
        
        var.append(Dis)
        inter_dis.append([tf, x, Dis])
    
    score = sum(var) / len(var)
    return [tf, round(score, 2), inter_dis, len(inter_dis)]



class FunctionTwo:

    def copula(path):
        print('START task functionTow.py')
        t1_start = process_time()

        sample_1_loom = os.path.join(path, 'sample_1_orgial-expression_matrix.loom')
        sample_2_loom = os.path.join(path, 'sample_2_orgial-expression_matrix.loom')

        lf = lp.connect(sample_1_loom, mode='r', validate=False)
        exprMat_sample_1 = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T

        lf = lp.connect(sample_2_loom, mode='r', validate=False)
        exprMat_sample_2 = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T

        modules = pd.read_csv(os.path.join(path, 'modules_all.csv'))

        # TODO pair by pair
        # TODO How size affect this
        copula_resl = []

        datalist = []
        for (columnName, columnData) in modules.iteritems():
            bool_series = pd.notnull(modules[columnName])
            datalist.append([list(modules[columnName][bool_series]), exprMat_sample_1, exprMat_sample_2])

        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(copulas, datalist)
            for r in results:
                copula_resl.append(r)
        end = time.perf_counter()
        print(f'Finished in {round(end - start, 2)} s')

        copula_s = pd.DataFrame(copula_resl)
        copula_s = copula_s.rename(columns={0: "TFs", 1: "Copula-score", 2: "inter_dis",3:"inter_Num"})

        print('********* here inter_dis check **********')
        print(copula_s)
        #print(copula_s.shape)
        #print(copula_s)
        inter_dis_list = copula_s["inter_dis"].tolist()
        # print('list')
        # print(inter_dis_list)
        inter_dis_path = os.path.join(path,'inter_dis.txt')
        with open(inter_dis_path, "wb") as fp:
            pickle.dump(inter_dis_list, fp)

        copula_s = copula_s.drop(columns=['inter_dis'], axis=1)

        # copula_s = copula_s.sort_values(by="Copula-score", ascending=False)
        copula_s["Rank-copula"] = copula_s["Copula-score"].rank(method='dense').astype(int)

        DC_scores = pd.read_csv(os.path.join(path, 'DC_scores.csv'))
        # DC_scores = DC_scores.sort_values(by="DC-score", ascending=False)
        DC_scores["Rank-DiffCo"] = DC_scores["DC-score"].rank(method='dense').astype(int)

        a = pd.merge(DC_scores, copula_s, on='TFs')
        # a["Score"] = a[['Rank-copula', 'Rank-DiffCo']].max(axis=1)
        # a['Score'] = (a['Rank-copula'] + a['Rank-DiffCo']) / 2

        x = a[['DC-score', 'Copula-score']]
        min_max_scaler = preprocessing.MinMaxScaler()
        x_scaled = min_max_scaler.fit_transform(x)
        dataset = pd.DataFrame(x_scaled, columns=['DC-score', 'Copula-score'])
        a['DC-normalized'] = dataset['DC-score']
        a['copula-normalized'] = dataset['Copula-score']
        # a['Score'] = (a['DC-normalized'] + a['copula-normalized'])
        # a['Score'] = a['copula-normalized']
        # a['Rank'] = a['Score'].rank(ascending=False)
        #a['Rank'] = a["Copula-score"].rank(method='dense').astype(int)
        a['Rank'] = a["DC-score"].rank(method='dense',ascending=False).astype(int)
        a.sort_values("Rank", inplace=True)
        print("LOOK INTO RANK: i did :)")
        print(a)

        a.to_csv(os.path.join(path, 'ranked_modules.csv'), index=None)
        t1_stop = process_time()
        print('time task functionTow.py:', t1_stop - t1_start)
        print('DONE task functionTow.py \n\n\n')
