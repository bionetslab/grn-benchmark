import pickle
from pyscenic.utils import modules_from_adjacencies
from distributed import LocalCluster, Client
import os
import glob
import argparse
import re
#from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.prune import prune2df, df2regulons
from pyscenic.utils import modules_from_adjacencies
import csv
import pandas as pd
import os
from os import chdir, getcwd
from os import path
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import sys
import subprocess
from time import process_time
from contextlib import contextmanager
import sys, os
from dask.distributed import LocalCluster, Client
from dask.diagnostics import ProgressBar
from pyscenic.utils import modules_from_adjacencies


# TODO all imports
# TODO give the user the choice of algorithm to use for GRN inference
# TODO revise the pysenic stuff
# TODO look into cisTarget in details
# TODO add dask to all the computing stuff (will make it fast) dask distribution


class scGRN:

    # TODO final number of workers
    def GRN_function(path_to_TFs, path_to_data, path_output_file_csv, path_to_matrix_csv, path_to_modules_p,
                     path_to_feather, f_motif_path, p_motifs_csv, p_regul_p):
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('START task scGRN.py: GRN_function')
        t1_start = process_time()
        workers_number = 40
        """This function take as input single cell expression matrices and return regulons"""
        # TODO we keep grnboost2 (a fast implementation og genie3 ) #options -m genie3 grnboost2
        # TODO down-regulated genes >>>> --all_modules 'yes'

        script_path = os.getcwd() + '/test_interface/scripts/see.py'
        subprocess.run(['python', script_path,
                        str(path_to_data),
                        str(path_to_TFs),
                        '--method', 'grnboost2',
                        '--output', str(path_output_file_csv),
                        '--num_workers', str(workers_number),
                        '--seed', '777'])

        adjacencies = pd.read_csv(path_output_file_csv, index_col=False, sep=',')

        f_final_loom = path_to_data
        lf = lp.connect(f_final_loom, mode='r', validate=False)
        exprMat = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T
        print('matrix size:' + str(exprMat.shape))
        exprMat.to_csv(path_to_matrix_csv, index=False)
        modules = list(modules_from_adjacencies(adjacencies, exprMat))
        print("***************************************")
        print('number of modules:' + str(len(modules)))
        print("***************************************")
        MODULES_FNAME = (path_to_modules_p)
        # saving
        with open(MODULES_FNAME, 'wb') as f:
            pickle.dump(modules, f)

        # TODO check if this works in case of many feather files (it suppose to work :))
        db_fnames = glob.glob(path_to_feather)

        def name(fname):
            return os.path.splitext(os.path.basename(fname))[0]

        dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
        print(dbs)

        f_motif_path = f_motif_path
        print(f_motif_path)
        print("prune2df")
        df = prune2df(dbs, modules, f_motif_path, num_workers=workers_number)
        # TODO do we  need this
        # df.to_csv(p_motifs_csv)
        print('df2regulons')
        regulons = df2regulons(df)
        # saving
        with open(p_regul_p, 'wb') as f:
            pickle.dump(regulons, f)
        regulons = [r.rename(r.name.replace('(+)', ' (' + str(len(r)) + 'g)')) for r in regulons]
        print("***************************************")
        print('Number of regulons:' + str(len(regulons)))
        print("***************************************")
        t1_stop = process_time()
        print('time:', (t1_stop - t1_start)/60,'min')

        print('DONE task scGRN.py: GRN_function')
        return regulons


class PathsG:

    # TODO what to save and what to not
    def paths_generator(paths, data_folder_sample, path_to_TFs, path_to_feather, f_motif_path):
        print('START task scGRN.py : PathsG')

        res = []
        for x in paths:
            path_to_data = x
            path_output_file_csv = x.split('.')[0] + '_' + 'result.csv'
            path_to_matrix_csv = x.split('.')[0] + '_' + 'exprMat.csv'
            path_to_modules_p = x.split('.')[0] + '_' + 'modules.p'
            p_motifs_csv = x.split('.')[0] + '_' + 'motifs.csv'
            p_regul_p = x.split('.')[0] + '_' + 'regulons.p'
            res.append(scGRN.GRN_function(path_to_TFs, path_to_data, path_output_file_csv, path_to_matrix_csv,
                                          path_to_modules_p, path_to_feather, f_motif_path, p_motifs_csv, p_regul_p))
        print('DONE task scGRN.py : PathsG')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        return res
