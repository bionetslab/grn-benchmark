import numpy as np
import pandas as pd
import scanpy as sc
from django.contrib.auth.models import User
from users.models import UserFiles
import loompy as lp
import os
from os import chdir, getcwd
from os import path
import random
from test_interface.tasks.scGRN import scGRN
from test_interface.tasks.scGRN import PathsG


class ComputeGrn:

    def get_matrices(userrr, n_sub, s_size, TFs_type):
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('START task compute_GRN.py')


        current_user = UserFiles.objects.get(user=userrr)
        adata = sc.read_h5ad(current_user.main_anndata)
        print("here")
        print('selection type')
        print(current_user.selectionType['selectionType'])
        print(adata)



        # TODO (not in this version)
        # work with highly variable genes or not
        # adata = adata[:, adata.var.highly_variable]

        if current_user.selectionType['selectionType'] == 'Lasso':
            cord = adata.obsm[current_user.map]
            cordn = []
            for x in cord:
                cordn.append([float('%.5f' % float(x[0])), float('%.5f' % float(x[1]))])
            xx1 = [round(float(w), 5) for w in current_user.sample1x_lasso]
            yy1 = [round(float(w), 5) for w in current_user.sample1y_lasso]
            xx2 = [round(float(w), 5) for w in current_user.sample2x_lasso]
            yy2 = [round(float(w), 5) for w in current_user.sample2y_lasso]
            include_cells_s1 = []
            for i in range(len(cordn)):
                for j in range(len(xx1)):
                    if cordn[i][0] == xx1[j] and cordn[i][1] == yy1[j]:
                        include_cells_s1.append(i)
            include_cells_s2 = []
            for i in range(len(cordn)):
                for j in range(len(xx2)):
                    if cordn[i][0] == xx2[j] and cordn[i][1] == yy2[j]:
                        include_cells_s2.append(i)

            sample_1 = adata[include_cells_s1]
            sample_2 = adata[include_cells_s2]

            print(sample_1)
            print(sample_2)

        else:
            sample1 = current_user.samples['sample1']
            sample2 = current_user.samples['sample2']
            anndatalist = []
            sample_1 = 0
            for x in sample1:
                anndatalist.append(adata[adata.obs[current_user.obs] == x, :])
            if len(anndatalist) > 1:
                sample_1 = anndatalist[0].concatenate(v for v in anndatalist[1:])
            else:
                sample_1 = adata[adata.obs[current_user.obs] == sample1[0], :]

            anndatalisttwo = []
            sample_2 = 0
            for y in sample2:
                anndatalisttwo.append(adata[adata.obs[current_user.obs] == y, :])
            if len(anndatalisttwo) > 1:
                sample_2 = anndatalisttwo[0].concatenate(v for v in anndatalisttwo[1:])
            else:
                sample_2 = adata[adata.obs[current_user.obs] == sample2[0], :]

            print(sample_1)
            print(sample_2)

        # data folder
        data_folder_all = os.getcwd() + '/media/' + str(userrr)

        def save_data(A, f_loom_path_scenic):
            # create basic row and column attributes for the loom file:
            row_attrs = {
                "Gene": np.array(A.var_names),
            }
            col_attrs = {
                "CellID": np.array(A.obs_names),
                "nGene": np.array(np.sum(A.X.transpose() > 0, axis=0)).flatten(),
                "nUMI": np.array(np.sum(A.X.transpose(), axis=0)).flatten(),
            }
            lp.create(f_loom_path_scenic, A.X.transpose(), row_attrs, col_attrs)

        # Save the main samples matrices files for copula analysis

        save_data(sample_1, path.join(data_folder_all, 'sample_1_orgial-expression_matrix.loom'))
        save_data(sample_2, path.join(data_folder_all, 'sample_2_orgial-expression_matrix.loom'))

        # TODO check is we need smoothing function

        def sub_sampling(sample, n_sub, s_size):
            """ This function do the subsampling of the inputted data matrices"""
            number_of_samples = int(n_sub)
            #sampling_size = int(s_size)
            sub_samples = []
            cells = list(sample.obs_names)
            print('original number of cells: '+str(len(cells)))
            sampling_size = int((int(s_size)*len(cells))/100)
            print('sampling_size_percentage : '+ str(sampling_size))

            # TODO add errors all around when needed if len(cells) < sampling_size: raise ValueError( "Invalid
            #  subsampling parameters. The data size {} is smaller than the sampling window  {}.".format(len(cells),
            #  sampling_size))

            # TODO look at the random pickup method it may make difference

            for i in range(number_of_samples):
                if len(cells)>sampling_size:
                    sub_samples.append(sample[sample.obs_names.isin(random.sample(cells, k=sampling_size))])
                else:
                    sub_samples.append(sample)
            return sub_samples

        # TODO add subsampling parameters selection to the frontend (DONE)

        sample_1_sub = sub_sampling(sample_1, n_sub, s_size)
        sample_2_sub = sub_sampling(sample_2, n_sub, s_size)
        print("*************sub-samples************")
        for x in sample_1_sub:
            print(x.shape)
        for x in sample_2_sub:
            print(x.shape)

        def auto_saving(samples, data_folder):
            if os.path.exists(data_folder) is False:
                os.mkdir(data_folder)
            paths = []
            index = 0
            for x in samples:
                f_loom_path = data_folder + 'sub_sample_' + str(index) + ".loom"
                paths.append(f_loom_path)
                save_data(x, f_loom_path)
                index = index + 1
            return paths

        data_folder_sample_1 = path.join(data_folder_all, 'data_folder_sample_1/')
        data_folder_sample_2 = path.join(data_folder_all, 'data_folder_sample_2/')

        sample_1_data_holder = auto_saving(sample_1_sub, data_folder_sample_1)
        sample_2_data_holder = auto_saving(sample_2_sub, data_folder_sample_2)

        # TODO add spices human/mouse to the parameters selection
        # TODO Download all dataset
        # folder location: static GRN_files
        # TFs_pl1573_yzr_3294_random.txt

        if current_user.main_TFs == '' and TFs_type == 'Human':
            print('Messages: You did not upload a TFs file (using all TFS). You are using Human TFs')
            path_to_TFs = os.getcwd() + "/test_interface/static/GRN_files/human/allTFs_hg38.txt"

        elif current_user.main_TFs == '' and TFs_type == 'Mouse':
            print('Messages: You did not upload a TFs file (using all TFS). You are using Mouse TFs')
            path_to_TFs = os.getcwd() + "/test_interface/static/GRN_files/Mouse/allTFs_mm.txt"

        else:
            print('Messages: You are using your own TFs file')
            path_to_TFs = current_user.main_TFs.path

        if os.stat(path_to_TFs).st_size == 0:
            raise ValueError(
                "Empty transcription factors file: please upload a file containing TFs in line text format.")

        # TODO download the data
        if TFs_type == 'Human':
            path_to_feather = os.getcwd() + '/test_interface/static/GRN_files/human/*mc9nr*.feather'
            f_motif_path = os.getcwd() + "/test_interface/static/GRN_files/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
        elif TFs_type == 'Mouse':
            path_to_feather = os.getcwd() + '/test_interface/static/GRN_files/Mouse/*mc9nr.feather'
            f_motif_path = os.getcwd() + "/test_interface/static/GRN_files/Mouse/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

        PathsG.paths_generator(sample_1_data_holder, data_folder_sample_1, path_to_TFs, path_to_feather, f_motif_path)
        PathsG.paths_generator(sample_2_data_holder, data_folder_sample_1, path_to_TFs, path_to_feather, f_motif_path)
        print('DONE task compute_GRN.py')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

