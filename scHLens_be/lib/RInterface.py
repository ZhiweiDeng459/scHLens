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

    r_script = open('./lib/function.R').read()
    robjects.r(r_script)
    devianceFS = robjects.globalenv['devianceFS']
    result = list(devianceFS(csc_matrix(X),var_index.to_list(),obs_index.to_list()))
    result = result[0:topGenes]
    return result

def SCTransform(X,var_index,obs_index,topGenes):
    '''
    X:adata.X
    var_index:adata.var.index
    obs_index:adata.obx.index
    topGenes: the choose top genes
    '''
    '''
    return : a list of chosen genes
    ''' 
    anndata2ri.activate()

    r_script = open('./lib/function.R').read()
    robjects.r(r_script)
    glmGamPoi = robjects.globalenv['MySCTransform']
    result = glmGamPoi(csc_matrix(X),var_index.to_list(),obs_index.to_list(),topGenes)
    
    # robjects.r['source']('./lib/function.R')
    # glmGamPoi = robjects.r['glmGamPoi']
    # result = glmGamPoi(csc_matrix(X),var_index.to_list(),obs_index.to_list())
    
    #转置
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

def CellChat(matrix,cellName,geneName,label,DatabaseType):
    '''
    matrix: A matrix of normlized data
    cellName: A list of cell name
    geneName: A list of gene name
    label: A list of cell's annonation
    DatabaseType: human or mouse
    '''
    '''
    return: A dict consist of the net count and net weight
    '''
    anndata2ri.activate()

    r_script = open('./lib/function.R').read()
    robjects.r(r_script)
    CellChat = robjects.globalenv['CellChat']
    if type(matrix) == np.ndarray:
        tempResult = CellChat(matrix,cellName,geneName,label,DatabaseType)
    elif type(matrix) == csr_matrix:
        tempResult = CellChat(csc_matrix(matrix),cellName,geneName,label,DatabaseType)
    elif type(matrix) == csc_matrix:
        tempResult = CellChat(matrix,cellName,geneName,label,DatabaseType)
    else:
        tempResult = None
    

    ## pack the result
    CellChatResult = {}
    types = ['count','weight']
    for t in types:
        clusters = list(tempResult.rx2(t).rx2('clusterList'))
        data = tempResult.rx2(t).rx2('data').tolist()
        dir_ = {}
        for i in range(0,len(clusters)):
            dir_[clusters[i]] = {}
            for j in range(0,len(clusters)):
                dir_[clusters[i]][clusters[j]] = data[i][j]
        CellChatResult[t] = dir_
        
    return CellChatResult

