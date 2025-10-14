'''
负责R函数的接口
'''
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import scanpy as sc
import anndata2ri
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
import rpy2.rinterface as rinterface
import numpy as np
from lib.utils import *
import subprocess
import json
import os

def scry(X,var_index,obs_index,topGenes):
    '''
    X:adata.X (numpy.ndarray)
    var_index:adata.var.index
    obs_index:adata.obx.index
    topGenes: the choose top genes
    '''
    '''
    return : a list of chosen genes
    ''' 
    anndata2ri.activate()

    ## remove numpy type
    # if isinstance(np.ndarray):
    X = X.astype(np.float64)        

    r_script = open('./lib/function.R', encoding="utf-8").read()
    robjects.r(r_script)
    devianceFS = robjects.globalenv['devianceFS']
    result = list(devianceFS(csc_matrix(X),var_index.to_list(),obs_index.to_list()))
    result = result[0:topGenes]
    return result

def SCTransform(X,var_index,obs_index,topGenes,JobId):
    '''
    X:adata.X
    var_index:adata.var.index
    obs_index:adata.obx.index
    topGenes: the choose top genes
    '''
    '''
    return : a list of chosen genes
    ''' 
    root_path = './job_attachment/' + JobId
    # omit_data
    save_matrix_for_R(get_dense_adata_X(X), obs_index, var_index, save_path=root_path+'/scTransform')
    params = {
        'topGenes':topGenes
    }
    with open(root_path + '/scTransform/params.json', "w", encoding="utf-8") as f:
        json.dump(params, f, ensure_ascii=False, indent=4)

    
    script_path = './lib/MySCtransform.R'
    subprocess.run(["Rscript", script_path] + [JobId], cwd=os.getcwd(),check=True)
    
    result = pd.read_csv(root_path + '/scTransform/result.csv', index_col=0)
    result = result.T
    
    return result

def slingshot(rd,cl):
    '''
    rd:rd:the DR matrix
    cl:A list of cell's annonation
    '''
    '''
    return : a dict consist of the lineages and the curves
    ''' 
    anndata2ri.activate()
    
    r_script = open('./lib/function.R').read()
    robjects.r(r_script)
    Slingshot = robjects.globalenv['Slingshot']
    tempResult = Slingshot(rd,cl)
    SlingshotResult = {}
    SlingshotResult['shape'] = {}
    for name,data in tempResult.rx2('lineages').items():
        SlingshotResult['shape'][name] = {'lineages':None,'curves':None}
        SlingshotResult['shape'][name]['lineages'] = list(data)
    for name,data in tempResult.rx2('curves').items():
        SlingshotResult['shape'][name]['curves'] = data.rx2('s').tolist()
    SlingshotResult['PseudotimeColor'] = list(tempResult.rx2('PseudotimeColor'))
    return SlingshotResult

def scCCESS_Kmeans(X,var_index,obs_index):
    '''
    X:adata.X
    var_index:adata.var.index
    obs_index:adata.obx.index
    '''
    '''
    return : the recommend clusters num
    '''
    
    anndata2ri.activate()

    r_script = open('./lib/function.R', encoding="utf-8").read()
    robjects.r(r_script)
    scCCESS_Kmeans_r = robjects.globalenv['scCCESS_Kmeans']
    result = scCCESS_Kmeans_r(csc_matrix(X),var_index.to_list(),obs_index.to_list())
    
    return result


# def CellChat(matrix,cellName,geneName,label,organism):
#     '''
#     matrix: A matrix of normlized data
#     cellName: A list of cell name
#     geneName: A list of gene name
#     label: A list of cell's annonation
#     organism: human or mouse
    
#     '''
#     '''
#     return: A dict consist of the net count and net weight
#     '''
#     anndata2ri.activate()

#     r_script = open('./lib/function.R', encoding="utf-8").read()
#     robjects.r(r_script)
#     CellChat = robjects.globalenv['CellChat']
#     if type(matrix) == np.ndarray:
#         tempResult = CellChat(matrix,cellName,geneName,label,organism)
#     elif type(matrix) == csr_matrix:
#         tempResult = CellChat(csc_matrix(matrix),cellName,geneName,label,organism)
#     elif type(matrix) == csc_matrix:
#         tempResult = CellChat(matrix,cellName,geneName,label,organism)
#     else:
#         tempResult = None
    

#     ## pack the result
#     CellChatResult = {}
    
#     ### 通讯基本信息
#     types = ['count','weight']
#     for t in types:
#         clusters = list(tempResult.rx2(t).rx2('clusterList'))
#         data = tempResult.rx2(t).rx2('data').tolist()
#         dir_ = {}
#         for i in range(0,len(clusters)):
#             dir_[clusters[i]] = {}
#             for j in range(0,len(clusters)):
#                 dir_[clusters[i]][clusters[j]] = data[i][j]
#         CellChatResult[t] = dir_
#     ### Interaction信息
#     types = ["source","target","ligand","receptor","prob"]
#     for t in types:
#         values = list(tempResult.rx2(t))
#         CellChatResult[t] = values
        
#     return CellChatResult


def CellChat(X,var_index,obs_index,label,organism,JobId):
    '''
    matrix: A matrix of normlized data
    var_index:adata.var.index
    obs_index:adata.obx.index
    label: A list of cell's annonation
    organism: human or mouse
    JobId: job id
    '''
    root_path = './job_attachment/' + JobId
    # omit_data
    save_matrix_for_R(get_dense_adata_X(X), obs_index, var_index, save_path=root_path+'/CellChat')
    params = {
        'label':label,
        'organism':organism
    }
    with open(root_path + '/CellChat/params.json', "w", encoding="utf-8") as f:
        json.dump(params, f, ensure_ascii=False, indent=4)

    
    script_path = './lib/CellChat.R'
    subprocess.run(["Rscript", script_path] + [JobId], cwd=os.getcwd(),check=True)
    
    # result = pd.read_csv(root_path + '/CellChat/result.csv', index_col=0)
    # result = result.T
    
    with open(root_path + '/CellChat/result.json', 'r', encoding='utf-8') as f:
        result = json.load(f)

    CellChatResult = {}
    ### 打包通讯基本信息
    types = ['count','weight']
    for t in types:
        CellChatResult[t] = {}
        clusters = result[t]['clusterList']
        for i in range(0,len(clusters)):
            CellChatResult[t][clusters[i]] = {}
            for j in range(0,len(clusters)):
                CellChatResult[t][clusters[i]][clusters[j]] = result[t]['data'][i][j]
    ### Interaction信息
    types = ["source","target","ligand","receptor","prob"]
    for t in types:
        CellChatResult[t] = result[t]

    return CellChatResult


    
    

def multik(X,var_index,obs_index,JobId):
    
    '''
    X:adata.X
    var_index:adata.var.index
    obs_index:adata.obx.index
    '''
    '''
    return : best k clusters
    ''' 
    
    root_path = './job_attachment/' + JobId
    # omit_data
    save_matrix_for_R(get_dense_adata_X(X), obs_index, var_index, save_path=root_path+'/multik')
    script_path = './lib/multik.R'
    subprocess.run(["Rscript", script_path] + [JobId], cwd=os.getcwd(),check=True)
    
    with open(root_path + "/multik/result.json", "r", encoding="utf-8") as f:
        data = json.load(f)
        best_k = data['best_k'][0]
    
    return best_k
