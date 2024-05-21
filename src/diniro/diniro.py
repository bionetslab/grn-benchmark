from django.shortcuts import render
from django.http import HttpResponse
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from users import views as user_viwes
from django.contrib.auth.models import User
from django.db import models
from users.models import UserFiles
from users.forms import UserFilesForm
from django.shortcuts import redirect
from wsgiref.util import FileWrapper
from test_interface.tasks.save_data import SaveData
from test_interface.tasks.compute_GRN import ComputeGrn
from test_interface.tasks.common import Common
from test_interface.tasks.scGRN import scGRN
from test_interface.tasks.scGRN import PathsG
from test_interface.tasks.function_one import FunctionOne
from test_interface.tasks.function_two import FunctionTwo
from test_interface.tasks.Filtering_processing import FilteringProcessing
from test_interface.tasks.go_tables import GoTerms
from test_interface.tasks.p_value import Pvalues
from test_interface.tasks.plotting import Plotting
from test_interface.tasks.gene_edges_distribution import GeneEdgesDist
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_protect
from django.core.files import File
import simplejson
import mimetypes
import pickle
import json
import numpy as np
import pandas as pd
import scanpy as sc
import json
import random
import os
import time
from os import path
import shutil


def rredirect(to, *args, permanent=False, **kwargs):
    redirect_response = redirect(to, *args, permanent=permanent, **kwargs)
    redirect_response['Location'] = '' + redirect_response['Location']
    return redirect_response

def home(request):
    print('##########################################################')
    print('********* view/home *********')
    context = {}
    path = os.getcwd() + '/media'
    now = time.time()

    for filename in os.listdir(path):
        filestamp = os.stat(os.path.join(path, filename)).st_mtime
        filecompare = now - 10000 * 3600
        if filestamp < filecompare and len(filename)==10:
            print(filename)
            print('files to delet')
            print(filename)
            pathremove = os.path.join(path, filename)
            try:
                shutil.rmtree(pathremove)
            except OSError as e:
                print("Error: %s : %s" % (pathremove, e.strerror))
            print('\n')

    print('##########################################################')
    return render(request, 'test_interface/home.html', context)


@csrf_protect
def upload(request):
    print(request.user)
    print('##########################################################')
    print('********* views.upload *********')
    if request.method == 'POST':
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        # creat user first
        user_viwes.creat_user(request)
        UserFiles.objects.create(user=request.user.username)
        model, created = UserFiles.objects.get_or_create(user=request.user.username)
        # user is created and logged in

        form = UserFilesForm(request.POST, request.FILES, instance=model)
        #print(request.FILES)
        isitdataexample = request.POST.get('thisdataexample')
        print(request)
        print(isitdataexample)
        example_data_is = 0
        if isitdataexample is not None:
            if int(isitdataexample) == 10 or int(isitdataexample)==20:
                #print(isitdataexample)
                example_data_is = 1

        if 1==1:

            if form.is_valid():
                print("valid")
                instance = form.save(commit=False)
                instance.save()
                example_data_is = 1
            if not form.is_valid():
                if isitdataexample is None:
                    print('not valid')
                    form = UserFilesForm()
                    context = {'form': form,
                               'FTV': 1
                               }
                    print('##########################################################1')
                    return render(request, 'test_interface/upload.html', context)
            if isitdataexample is not None and int(isitdataexample) == 10:
                #('&&&&&&&&&&&&&&&&&&')
                instance, created = UserFiles.objects.get_or_create(user=request.user.username)
                pathtofile = os.getcwd()+'/media/usethis350genes.h5ad'
                #print(pathtofile)
                f = open(pathtofile, "rb")
                instance.main_anndata = File(f)
                instance.sessionname = 'Satija et al.(2015)-data'
                instance.save()
            if isitdataexample is not None and int(isitdataexample) == 20:
                #('&&&&&&&&&&&&&&&&&&')
                instance, created = UserFiles.objects.get_or_create(user=request.user.username)
                pathtofile = os.getcwd()+'/media/test4.h5ad'
                #print(pathtofile)
                f = open(pathtofile, "rb")
                instance.main_anndata = File(f)
                instance.sessionname = 'Guo X et al.(2018)-data'
                instance.save()

            # Creat data for plotting
            userrr = request.user.username
            N_Cells, N_Genes, N_Groups, points, leb, cc, plotting_maps, observations = SaveData.read_data(userrr)
            unique_leb = list(set(leb))

            plotting_data = SaveData.legend_check(unique_leb, points, leb, cc)

            instance, created = UserFiles.objects.get_or_create(user=request.user.username)
            #print(N_Cells)
            print(instance.main_anndata)
            print(instance.sessionname)
            sessionnamegiven = instance.sessionname
            userr = request.user.username
            if instance.main_TFs =='':
                print('You are using all known TFs')
            context = {'form': form,
                       'isitdataexample': example_data_is,
                       'instance': instance,
                       'leb': leb,
                       'unique_leb': unique_leb,
                       'plotting_data': plotting_data,
                       'plotting_maps': plotting_maps,
                       'observations': observations,
                       'N_Cells': N_Cells,
                       'N_Genes': N_Genes,
                       'N_Groups': N_Groups,
                       'userr':userr,
                       'sessionname':sessionnamegiven
                       }
            print('##########################################################2')
            return render(request, 'test_interface/upload.html', context)


    else:
        form = UserFilesForm()
        context = {'form': form}
        print('##########################################################3')
        return render(request, 'test_interface/upload.html', context)


def selection1(request):
    print('##########################################################')
    print('********* view.selection I *********')
    print(request.user)

    map = request.POST.get('map')
    obs = request.POST.get('obs')
    userrr = request.POST.get('user')
    reset = request.POST.get('reset')
    print('reset :'+str(reset))

    print(map)
    print(obs)
    print(userrr)

    model, created = UserFiles.objects.get_or_create(user=userrr)
    print(model.main_anndata)

    if int(reset) == 1:
        print('here')

        cc_save = os.getcwd() + '/media/' + str(userrr) + '/plotting_data.txt'
        with open(cc_save, "rb") as fp:
            cc1 = pickle.load(fp)

        model.obs = obs
        model.map = map
        model.save()
        # Creat data for plotting
        N_Cells, N_Genes, N_Groups, points, leb, cc, plotting_maps, observations = SaveData.read_data(userrr, map, obs)
        unique_leb = list(set(leb))

        plotting_data = SaveData.legend_check(unique_leb, points, leb, cc1)

        context = {'leb': leb,
                   'unique_leb': unique_leb,
                   'plotting_data': plotting_data,
                   'plotting_maps': plotting_maps,
                   'observations': observations
                   }
        print('##########################################################')
        return JsonResponse(context)


    if model.obs == obs and model.map == map and int(reset) == 0:
        print('skip')
        # model.obs = obs
        # model.map = map
        # model.save()
        # # Creat data for plotting
        # N_Cells, N_Genes, N_Groups, points, leb, cc, plotting_maps, observations = SaveData.read_data(userrr, map, obs)
        # unique_leb = list(set(leb))
        # cc =
        # plotting_data = SaveData.legend_check(unique_leb, points, leb, cc)

        context = {
                   # 'leb': leb,
                   # 'unique_leb': unique_leb,
                   # 'plotting_data': plotting_data,
                   # 'plotting_maps': plotting_maps,
                   # 'observations': observations
                   }
        print('##########################################################')
        return JsonResponse(context)

    else:
        print('dont skip')
        model.obs = obs
        model.map = map
        model.save()
        # Creat data for plotting
        N_Cells, N_Genes, N_Groups, points, leb, cc, plotting_maps, observations = SaveData.read_data(userrr, map, obs)
        unique_leb = list(set(leb))
        plotting_data = SaveData.legend_check(unique_leb, points, leb, cc)

        cc_save = os.getcwd() + '/media/' + str(userrr) + '/plotting_data.txt'
        with open(cc_save, "wb") as fp:
            pickle.dump(cc, fp)

        context = {'leb': leb,
                   'unique_leb': unique_leb,
                   'plotting_data': plotting_data,
                   'plotting_maps': plotting_maps,
                   'observations': observations
                   }
        print('##########################################################')
        return JsonResponse(context)


def selection2(request):
    print('##########################################################')
    print('********* view.selection II *********')
    print(request.GET.get)

    sample1 = request.GET.get('sample_1')
    sample2 = request.GET.get('sample_2')
    userr =   request.GET.get('user')


    Sample1Nmae = request.GET.get('s1name')
    Sample2Nmae = request.GET.get('s2name')
    model, created = UserFiles.objects.get_or_create(user=userr)
    model.sample1Name = Sample1Nmae
    model.sample2Name = Sample2Nmae
    model.save()
    print("************* Samples ID **********************")
    print(Sample1Nmae)
    print(Sample2Nmae)

    sample1 = sample1.split(',')
    sample2 = sample2.split(',')

    print(sample1)
    print(sample2)
    #print((list(sample1[0])))
    if len(list(sample1[0])) != 0 or len(list(sample1[0]))!=0:

        model, created = UserFiles.objects.get_or_create(user=userr)
        model.samples = {'sample1': sample1,
                         'sample2': sample2,
                         }
        model.selectionType = {'selectionType': 'Clusters'}
        model.save()

    context = {}
    print('##########################################################')
    return render(request, 'test_interface/upload.html')

@csrf_protect
def lasso_selecton(request):
    print('##########################################################')
    userrr = request.POST.get('user')
    currentsample = request.POST.get('drawsam')
    print(userrr)
    print(currentsample)
    print('********* Lasso.Selection *********')

    sample1 = 0
    sample2 = 0

    print('lasso_selecton')
    x = request.POST
    xa = list(dict(x).values())[1]
    ya = list(dict(x).values())[2]


    l = str(len(xa))
    print('Number of selected cells')
    print(len(xa))
    #print(xa)
    print(len(ya))
    #print(ya)
    print('End')

    model, created = UserFiles.objects.get_or_create(user=userrr)
    # print(len(model.sample1x_lasso))
    if int(currentsample)==1:
        print('adding to sample I')
        model.sample1x_lasso = xa
        model.sample1y_lasso = ya
        model.selectionType = {'selectionType': 'Lasso'}
        model.save()
        sample1 = 1
        sample2 = 0


    model, created = UserFiles.objects.get_or_create(user=userrr)
    # print(len(model.sample2x_lasso))
    if int(currentsample)==2:
        model.sample2x_lasso = xa
        model.sample2y_lasso = ya
        model.save()
        sample2 = 1
        sample1 = 0


    model, created = UserFiles.objects.get_or_create(user=userrr)
    print(len(model.sample1x_lasso))
    print(len(model.sample2x_lasso))

    #print(model.sample1x_lasso)
    #print(model.sample2x_lasso)

    context = {'lengthx': l,
               'sample1': sample1,
               'sample2': sample2}

    print('##########################################################')
    return JsonResponse(context)


@csrf_protect
def lasso_reset(request):
    print('##########################################################')
    print('********* Lasso.Reset *********')
    userrr = request.POST.get('user')

    default = {'foo': 'bar'}
    model, created = UserFiles.objects.get_or_create(user=userrr)
    model.sample1x_lasso = default
    model.sample1y_lasso = default
    model.sample2x_lasso = default
    model.sample2y_lasso = default
    model.save()

    context = {'reset': 1}
    print('##########################################################')
    return render(request, 'test_interface/upload.html',context)

    #return JsonResponse(context)


def result(request):
    print('##########################################################')
    print('********* view.result *********')

    if request.method == 'POST':
        n_sub = request.POST.get('n_sub')
        s_size = request.POST.get('s_size')
        TFs_type = request.POST.get('Species_selection')
        occurrence_threshold = request.POST.get('Occurrence_threshold')
        size_threshold_max = request.POST.get('Size_threshold_max')
        size_threshold_min = request.POST.get('Size_threshold_min')
        number_of_modules = request.POST.get('number_modules')
        pval_cutoff = request.POST.get('pval')
        pval_cutoff = float(pval_cutoff)
        userrr = request.POST.get('user')
        print(userrr)
        # TODO add this to parameters
        print('************ Paramters ***************')
        print(n_sub)
        print(s_size)
        print(TFs_type)
        print(size_threshold_max)
        print(size_threshold_min)
        print(occurrence_threshold)
        print(number_of_modules)
        print(pval_cutoff)
        # TODO remove this from some functions
        size_thresholds = 0

        # data_processing_and_GRN.ipynb
        # TODO allow this when you are done
        ComputeGrn.get_matrices(userrr, n_sub, s_size, TFs_type)
        # only saves tha data to the user folder

        # reg_preprocessing_and_filtering.ipynb
        #current_user = request.user.username
        current_user = userrr

        path_sample_1 = os.getcwd() + '/media/' + str(current_user) + '/data_folder_sample_1/'
        path_sample_2 = os.getcwd() + '/media/' + str(current_user) + '/data_folder_sample_2/'
        # TODO add parameters to the front end

        sample_1 = FilteringProcessing.merging_and_filtring(path_sample_1, occurrence_threshold, size_thresholds)
        sample_2 = FilteringProcessing.merging_and_filtring(path_sample_2, occurrence_threshold, size_thresholds)
        data_folder_all = os.getcwd() + '/media/' + str(current_user)
        #print('path' + str(data_folder_all))
        sample_1.to_csv(path.join(data_folder_all, 'sample_grn_s1.csv'), index=False)
        sample_2.to_csv(path.join(data_folder_all, 'sample_grn_s2.csv'), index=False)


        # function1.ipynb
        FunctionOne.functionOne(data_folder_all,size_thresholds)

        # function_2_copula.ipynb
        FunctionTwo.copula(data_folder_all)


        # p-values and cutoff filtering :

        p_values = Pvalues.p_values(data_folder_all)

        # after p value filtering

        GeneEdgesDist.update(data_folder_all,pval_cutoff)


        # plotting

        Nlist, Elist, nlist, elist, n_combined, e_combined = Plotting.plot_network(data_folder_all, number_of_modules,pval_cutoff,int(size_threshold_min),int(size_threshold_max))

        nodes = os.getcwd() + '/media/' + str(current_user) + '/nodes.txt'
        edges = os.getcwd() + '/media/' + str(current_user) + '/edges.txt'

        with open(nodes, "wb") as fp:
            pickle.dump(nlist, fp)

        with open(edges, "wb") as fp:
            pickle.dump(elist, fp)

        nodes2 = os.getcwd() + '/media/' + str(current_user) + '/nodes2.txt'
        edges2 = os.getcwd() + '/media/' + str(current_user) + '/edges2.txt'

        with open(nodes2, "wb") as fp:
            pickle.dump(Nlist, fp)

        with open(edges2, "wb") as fp:
            pickle.dump(Elist, fp)

        table = GoTerms.go_terms(data_folder_all,pval_cutoff,int(size_threshold_min),int(size_threshold_max))
        context = {}
        print('##########################################################')
        return render(request, 'test_interface/result.html', context)

    else:

        context = {}
        return render(request, 'test_interface/upload.html', context)



def result_page(request):

    #print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n")
    print('********* view/result page *********')
    #print(request)
    userrr = request.GET.get('user')
    #print(userrr)
    #print("POST")
    if request.method == 'POST':
        userrr = request.POST.get('userr')
    #print(userrr)
    #print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n")


    current_user = userrr
    data_folder_all = os.getcwd() + '/media/' + str(current_user)

    #current_user = 'oTsuRNfrnhrTlcqXmXbzcbf0lHAJdnAD'
    #print(current_user)

    instance, created = UserFiles.objects.get_or_create(user=current_user)
    sessionnamegiven = instance.sessionname

    figure_url = '/media/' + str(current_user) + '/common_interactions.png'
    figure1_url = '/media/' + str(current_user) + '/ranking_plot.png'

    go_file = os.getcwd() + '/media/' + str(current_user) + '/GO_terms.txt'
    if os.path.isfile(go_file):
        print("File exist")
        with open(go_file, "rb") as fp:
            table = pickle.load(fp)
    else:
        print("File not exist")
        table = []




    nodes = os.getcwd() + '/media/' + str(current_user) + '/nodes.txt'
    with open(nodes, "rb") as fp:
        n_combined = pickle.load(fp)

    edges = os.getcwd() + '/media/' + str(current_user) + '/edges.txt'
    with open(edges, "rb") as fp:
        e_combined = pickle.load(fp)

    nodes2 = os.getcwd() + '/media/' + str(current_user) + '/nodes2.txt'
    with open(nodes2, "rb") as fp:
        NODES = pickle.load(fp)

    edges2 = os.getcwd() + '/media/' + str(current_user) + '/edges2.txt'
    with open(edges2, "rb") as fp:
        EDGES = pickle.load(fp)

    #TODO    remove
    modules = []
    for x in table:
        modules.append(x[1])

    #TODO it start from here:
    #print('###########################################################################################################\n\n\n')
    SearchedGenes = request.POST.get('SearchedGenes')
    k = request.POST.get('k')
    view = request.POST.get('view')
    #print('###########################################################################################################\n\n\n')




    # if request.POST.get('SearchedGenes') == 'NoGene':
    #
    #     #print('This for the table networks')
    #     ID = request.POST.get('ModuleId')
    #     #print(ID)
    #     inc = []
    #
    #     for i in range(len(n_combined)):
    #           if len(n_combined[i]) !=0:
    #               if str(n_combined[i][0]['label']) == ID:
    #                 inc.append(i)
    #
    #     n = []
    #     e = []
    #     for i in inc:
    #         n += n_combined[i]
    #         e += e_combined[i]
    #     context = {'n_combined': n,'e_combined': e}
    #     return JsonResponse(context)
    Nn = []
    Ee = []
    for x in modules:
        inc = []
        for i in range(len(n_combined)):
            if len(n_combined[i]) != 0:
                if str(n_combined[i][0]['label']) == x:
                    inc.append(i)

        n = []
        e = []
        for i in inc:
            n += n_combined[i]
            e += e_combined[i]
        Nn.append(n)
        Ee.append(e)



    
    if k is not None:
        if int(k)==1:
            print('SearchedGenes and k and view')
            print(SearchedGenes)
            print(k)
            print(view)
            # Get the data:
            gene1Data, gene2Data, geneBothData, edges1Data, edges2Data, edgesBothData, As, At, Bs, Bt, Cs, Ct = GeneEdgesDist.vann_diagrams(data_folder_all)
            A = gene1Data
            B = gene2Data
            C = geneBothData
            A_edges = edges1Data
            B_edges = edges2Data
            C_edges = edgesBothData

            print("venn diagram changes ")
            GenesOfIntrest1 = SearchedGenes.split(',')
            GenesOfIntrest = []
            for z in GenesOfIntrest1:
                GenesOfIntrest.append(z.strip())

            print(GenesOfIntrest)
            print('view')
            print(view)

            AA = [i for i in GenesOfIntrest if i in A]
            BB = [i for i in GenesOfIntrest if i in B]
            CC = [i for i in GenesOfIntrest if i in C]

            #edges :
            AA_edges = []
            BB_edges = []
            CC_edges = []
            FoundGenes = []

            for g in GenesOfIntrest:
                for e in A_edges:
                    if g == e.split('>')[0] or g == e.split('>')[1]:
                        AA_edges.append(e.split('>')[0]+str('-->')+e.split('>')[1])
                        FoundGenes.append(g)

            for g in GenesOfIntrest:
                for e in B_edges:
                    if g == e.split('>')[0] or g == e.split('>')[1]:
                        BB_edges.append(e.split('>')[0]+str('-->')+e.split('>')[1])
                        FoundGenes.append(g)

            for g in GenesOfIntrest:
                for e in C_edges:
                    if g == e.split('>')[0] or g == e.split('>')[1]:
                        CC_edges.append(e.split('>')[0]+str('-->')+e.split('>')[1])
                        FoundGenes.append(g)

            print('search edges')
            print(AA_edges)
            print(BB_edges)
            print(CC_edges)
            NN_edges = [i for i in GenesOfIntrest if i not in FoundGenes]


            all = AA+BB+CC
            NN = [i for i in GenesOfIntrest if i not in all]

            #toggle:
            if int(view) == 1:
                AA = AA_edges
                BB = BB_edges
                CC = CC_edges
                NN = NN_edges

            print("%%%%%%%%")

            print(AA)
            print(BB)
            print(CC)
            print(NN)

            CCC = json.dumps(CC)
            AAA = json.dumps(AA)
            BBB = json.dumps(BB)
            NNN = json.dumps(NN)

            context = {'user':json.dumps(current_user),'la': len(AA),'lb': len(BB),'lc':len(CC),'ga': AAA,'gb': BBB,'gc': CCC,'gn': NNN,'ln': len(NN),'a': A,'b': B,'c': C,'As': As,'at':At,'bs': Bs,'bt': Bt,'cs': Cs,'ct': Ct,}

            return JsonResponse(context)

    elif request.method != 'POST':
        gene1Data, gene2Data, geneBothData, edges1Data, edges2Data, edgesBothData, As, At, Bs, Bt, Cs, Ct = GeneEdgesDist.vann_diagrams(data_folder_all)
        A = gene1Data
        B = gene2Data
        C = geneBothData
        A_edges = edges1Data
        B_edges = edges2Data
        C_edges = edgesBothData

        CCC = json.dumps(A_edges)
        AAA = json.dumps(B_edges)
        BBB = json.dumps(C_edges)


        print('No condition')
        userr = UserFiles.objects.get(user=current_user)
        s1name = userr.sample1Name
        s2name = userr.sample2Name
        print(s1name)
        print(s2name)

        n = n_combined[0]
        e = e_combined[0]

        #print(EDGES)
        print('size of nodes : '+str(len(NODES)))
        print('size of edges : '+str(len(EDGES)))
        NODES = NODES

        context = {'user': json.dumps(current_user),'table': table,'modules': json.dumps(modules),'n_combined': json.dumps(Nn),'e_combined': json.dumps(Ee), 'EDGES': EDGES,'NODES': NODES,'figure1_url': figure1_url, 'url': figure_url,'s1n': s1name,'s2n': s2name,
            'la': len(A),'lb': len(B),'lc': len(C),'ga': AAA,'gb': BBB,'gc': CCC,'a': A,'b': B,'c': C,'sessionname': sessionnamegiven,'As': As,'at': At,'bs': Bs,'bt': Bt,'cs': Cs,'ct': Ct}

        return render(request, 'test_interface/result.html', context)



def about(request):
    return render(request, 'test_interface/about.html')

# def home2(request):
#     return render(request, 'test_interface/home.html')

from django.utils.encoding import force_text, smart_str

def download(request):
    userr =   request.GET.get('user')
    file_path = os.getcwd() + '/media/' + userr + '/table.csv'
    print(file_path)
    file_wrapper = FileWrapper(open(file_path,'rb'))
    file_mimetype = mimetypes.guess_type(file_path)
    response = HttpResponse(file_wrapper, content_type=file_mimetype )
    response['X-Sendfile'] = file_path
    response['Content-Length'] = os.stat(file_path).st_size
    response['Content-Disposition'] = 'attachment; filename=%s' % smart_str('table.csv')
    return response

