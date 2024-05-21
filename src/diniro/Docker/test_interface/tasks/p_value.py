import anndata
import loompy as lp
import pandas as pd
import numpy as np
from copulae import EmpiricalCopula
from itertools import permutations
from copulae.core import pseudo_obs
import scipy
import random
import pickle
import os
import time
import math
from statsmodels.stats.multitest import multipletests

class Pvalues:


    def p_values(path):
        print('p_value.py')
        start = time.perf_counter()

        print("*********** p-value *****************************")
        sample_1_loom = os.path.join(path, 'sample_1_orgial-expression_matrix.loom')
        sample_2_loom = os.path.join(path, 'sample_2_orgial-expression_matrix.loom')

        lf = lp.connect(sample_1_loom, mode='r', validate=False)
        exprMat_sample_1 = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T

        lf = lp.connect(sample_2_loom, mode='r', validate=False)
        exprMat_sample_2 = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T


        def copula(df11, df22):
            cop1 = EmpiricalCopula(pseudo_obs(df11), ties='max')
            ec1 = cop1.cdf(pseudo_obs(df11))

            cop2 = EmpiricalCopula(pseudo_obs(df22), ties='max')
            ec2 = cop2.cdf(pseudo_obs(df22))

            Dis = 1-scipy.stats.ks_2samp(ec1, ec2)[0]

            return Dis

        C1_original = exprMat_sample_1.to_numpy()
        C2_original = exprMat_sample_2.to_numpy()

        print(C1_original.shape)
        print(C2_original.shape)

        c1_shape = C1_original.shape[0] * C1_original.shape[1]
        c2_shape = C2_original.shape[0] * C2_original.shape[1]

        C1 = C1_original.reshape(c1_shape, )
        C2 = C2_original.reshape(c2_shape, )

        print(C1.shape, C2.shape)

        np.random.shuffle(C1)
        np.random.shuffle(C2)

        C1 = C1.reshape(C1_original.shape)
        C2 = C2.reshape(C2_original.shape)

        print(C1.shape)
        print(C2.shape)

        C1_random = pd.DataFrame(data=C1)
        C2_random = pd.DataFrame(data=C2)

        if C1_random.shape[1] == C2_random.shape[1]:

            a = list(permutations(range(C2_random.shape[1]), 2))
            random.shuffle(a)
            #TODO change to 10000
            N = 10000
            random_pairs = a[:N]
            nulldis = []
            zeros_counts1 = []
            zeros_counts2 = []
            for x in random_pairs:
                df_c1 = C1_random[[x[0], x[1]]]
                df_c2 = C2_random[[x[0], x[1]]]

                z1 = list(df_c1.astype(bool).sum(axis=0))
                z2 = list(df_c2.astype(bool).sum(axis=0))

                zeros_counts1.append(df_c1.shape[0]-z1[0])
                zeros_counts1.append(df_c1.shape[0]-z1[1])
                zeros_counts2.append(df_c2.shape[0]-z2[0])
                zeros_counts2.append(df_c2.shape[0]-z2[1])


                #zeros1 = list(df_c1.isin([0]).sum())
                #zeros2 = list(df_c2.isin([0]).sum())
                #TODO filter this its important

                # for w in zeros1:
                #     if w/df_c1.shape[0] >= 0.3:
                #         print('Error too many zeros when randomising')
                # for y in zeros2:
                #     if y/df_c2.shape[0] >= 0.3:
                #         print('Error too many zeros when randomising')

                nulldis.append(copula(df_c1, df_c2))
            print('some zeros count after randomisation')
            # print(max(zeros_counts1))
            # for w in zeros_counts1:
            #     if w/C1.shape[0] >= 0.3:
                    #print('Error too many zeros when randomising')
            # for y in zeros_counts2:
            #     if y/C2.shape[0] >= 0.3:
                    #print('Error too many zeros when randomising')


        else:
            print('They should have same number of genes')

        with open(os.path.join(path, 'inter_dis.txt'), "rb") as fp:
            inter_dis = pickle.load(fp)

        inter_dis_all = []
        for x in inter_dis:
            for y in x:
                inter_dis_all.append(y)
        #print('***** inter_dis_all  ********************')
        #print(len(inter_dis_all))
        #print(inter_dis_all)
        p_value = []
        pval = []
        pval2 = []
        qval = []
        for true in inter_dis_all:
            s = sum(i >= true[2] for i in nulldis)
            s2 = sum(i <= true[2] for i in nulldis)
            p = s / N
            p2 = s2 / N
            #TODO just add it
            #p_enhanced = min(abs(1 - p), abs(0 - p))
            p_value.append([true[0], true[1], p])
            pval.append(p)
            pval2.append(p2)

        c = 0
        for x in pval:
            if x < 0.05:
                c = c + 1
        c1 = 0
        for x2 in pval2:
            if x2 < 0.05:
                c1 = c1 + 1

        print('this the used p value ' + str(c))
        print('this the other side of p value ' + str(c1))

        print("**************************** q-values ******************************")
        qvals = multipletests(pval, method="fdr_bh")[1]
        #print(qvals)
        print("**************************** p-values ******************************")
        #print(p_value[:5])
        p_vales_path = os.path.join(path,'p_values.txt')
        with open(p_vales_path, "wb") as fp:
            pickle.dump(p_value, fp)

        end = time.perf_counter()
        print(f'Finished in {round(end - start, 2)} s')
        print('\n\n\n')


