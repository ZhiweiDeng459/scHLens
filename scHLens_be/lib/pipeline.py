from lib.utils import *
import scanpy as sc
import numpy as np
import pandas as pd
import umap
import time
import random

import squidpy as sq
import os.path
import sc3s
# import scanorama
from lib.RInterface import *
from lib.corr import DIcorrect
import openTSNE
from scipy.sparse import csr_matrix,issparse
from scipy.spatial.distance import pdist, squareform
import cosg as cosg
import importlib
import gseapy as gp
import json
import seaborn as sns
from lib.vars import color_list
from flask_socketio import SocketIO, send,emit,join_room,leave_room
from sklearn.neighbors import NearestNeighbors
import copy
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

importlib.reload(cosg)


class pipelineException(Exception):
    def __init__(self,location='unKnown',advice='unKnown',message='unKnown'):
        self.location = location
        self.advice = advice
        self.message = message
'''

'''



'''
Pipeline
'''

def globalPipeline(params):
    
    ## JobId
    JobId = params['JobId']


    emit('get_pipeline_schedule',{'status':'Loading Dataset','percentage':0.2/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## Info: Loading Dataset

    ## 读文件
    try:
        adata,hasNan = readData(params['dataset'],JobId)
        if hasNan:
            emit('general_info',{
                'type':'warning',
                'title':'Clean invalid Values',
                'content':'The data set contains invalid values (such as null values), which have been automatically converted to 0.',
            },to=f'general_info_{JobId}',namespace='/')

    except Exception as e:
        raise pipelineException(
            location='Reading Dataset',
            advice='The dataset cannot be loaded, possibly because the dataset has issues (such as corruption or formatting errors), or because the required dataset files are missing.',
            message=str(e))

    ## ViewId
    try:
        ViewId = initView(JobId)
    except Exception as e:
        raise pipelineException(
            locaition='Initializing View',
            advice='Due to a server issue, the view initialization failed.',
            message=str(e))
    
    ## 填入参数
    adata.uns['params'] = params
    adata.uns['JobId'] = JobId
    adata.uns['ViewId'] = ViewId
    adata.uns['ParentId'] = params['ParentId']

    ## 采样
    emit('get_pipeline_schedule',{'status':'Sampling','percentage':0.5/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Sampling

    if 'SP' in adata.uns['params']:
        adata = SP(adata)

    emit('get_pipeline_schedule',{'status':'Quality Control','percentage':1.5/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Quality Control

    ## 质量控制 对每个数据集进行过滤
    if 'QC' in adata.uns['params']:
        adata = QC(adata)
    
    ## 数据融合
    if 'DI' in adata.uns['params'] and len(adata.uns['params']['DI']['Datasets']) != 0:
        adata = DI(adata)
    
    ## 如果不自带label，那么给予统一的默认label id：c_0
    if 'label' not in adata.obs:
        adata.obs['label'] = pd.Series(['c_0' for i in range(len(adata.obs))],dtype='category',index=adata.obs.index)
    
    localPointAdata = adata.copy()
    
    ## 存入原始矩阵 TODO 和数据融合配合
    adata.uns['count'] = adata.copy() #是size被QC改变之后的uns['count']

    emit('get_pipeline_schedule',{'status':'Transformation','percentage':2.0/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Transformation

    ## 变换 （FS发生在TS中）
    if 'TS' in adata.uns['params']:
        adata = TS1(adata)

    saveCache(adata,adata.uns['JobId'],adata.uns['ViewId'],'Query')

    ## 特征选择
    if 'FS' in adata.uns['params']:
        adata = FS(adata)

    if 'TS' in adata.uns['params']:
        adata = TS2(adata)

    emit('get_pipeline_schedule',{'status':'Neighboring','percentage':3.0/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Neighboring

    # ## neighbors
    # if 'NB' in adata.uns['params']:
    #     adata = NB(adata)

    emit('get_pipeline_schedule',{'status':'Visualization','percentage':4.0/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Visualization

    ## 降维
    if 'DR' in adata.uns['params']:
        adata = DR(adata)
    else:
        adata.obsm['embedding'] = np.zeros((len(adata.obs),2)) 

    emit('get_pipeline_schedule',{'status':'Clustering','percentage':5.0/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Visualization

    ## 聚类
    if 'CL' in adata.uns['params']:
        adata = CL(adata)
    # else:
    #     type_truth_key = 'cell_type1'
    #     adata.obs['label'] = adata.obs[type_truth_key] ##针对muraro数据集的特殊处理

    ## 保存聚类结果到localPoint
    ## 保存“供local pipeline进行划分时加载的缓存”
    localPointAdata.obs['label'] = adata.obs.label
    saveCache(localPointAdata, JobId, ViewId, 'localPoint')# uns为空，label只有
    

    # ## 过滤单样本的聚类，防止影响后续
    # adata = clearSimpleSizeCluster(adata)

    # ## 轨迹推断
    # if 'TI' in adata.uns['params']:
    #     adata = TI(adata)
    


    emit('get_pipeline_schedule',{'status':'Calculating DEGs','percentage':6.0/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Calculating Marker Genes


    ## 计算标志基因
    if 'MK' in adata.uns['params'] and len(adata.obs['label'].cat.categories) > 1:
        adata = MK(adata)
        adata.uns['init_raw_marker'] = adata.uns['raw_marker']


    emit('get_pipeline_schedule',{'status':'Cell Chat','percentage':6.5/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Assembling Results


    ## 细胞通讯
    if 'CC' in adata.uns['params']:
        adata = CC(adata)
        
    emit('get_pipeline_schedule',{'status':'Assembling Results','percentage':7.0/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Assembling Results

    ## 保存 “最终缓存”
    saveCache(adata,adata.uns['JobId'],adata.uns['ViewId'],'Last')

    ## 构建保存metaData
    metaData = generateMetaDataFromAdata(adata)
    saveViewMetaData(JobId, ViewId, metaData)

    ## 构建保存Tree
    saveToTree(adata)

    return adata


def localPipeline(params):
    ## JobId
    JobId = params['JobId']
    ## ParentId
    ParentId = params['ParentId']
    ## ViewId
    ViewId = initView(JobId)

    emit('get_pipeline_schedule',{'status':'Loading Cache','percentage':0.5 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Loading Cache

    ## 读取文件    
    adata = readCache(JobId, ParentId, 'localPoint')

    ## 局部化
    raw_celllist = adata.obs.index.tolist()
    chosenCells = sorted(params['type']['local']['chosenCells'],key=lambda x: raw_celllist.index(x)) #对过滤的细胞进行排序，防止造成随机性的问题
    adata = adata[chosenCells,:].copy()

    
    ## 填入参数
    adata.uns['params'] = params
    adata.uns['JobId'] = JobId
    adata.uns['ViewId'] = ViewId
    adata.uns['ParentId'] = ParentId


    emit('get_pipeline_schedule',{'status':'Quality Control','percentage':0.7/7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Quality Control

    ## 质量控制 对数据集进行过滤
    if 'QC' in adata.uns['params']:
        adata = QC(adata)


    ## 如果不自带label，那么给予统一的默认label id：c_0
    if 'label' not in adata.obs:
        adata.obs['label'] = pd.Series(['c_0' for i in range(len(adata.obs))],dtype='category',index=adata.obs.index)

    localPointAdata = adata.copy()

    ## 存入原始矩阵 TODO 和数据融合配合
    adata.uns['count'] = adata.copy()

    emit('get_pipeline_schedule',{'status':'Transformation','percentage':1.0 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Transformation


    ## 变换 （FS发生在TS中）
    if 'TS' in adata.uns['params']:
        adata = TS1(adata)

    saveCache(adata,adata.uns['JobId'],adata.uns['ViewId'],'Query')

    ## 特征选择
    if 'FS' in adata.uns['params']:
        adata = FS(adata)

    if 'TS' in adata.uns['params']:
        adata = TS2(adata)
        
        
    emit('get_pipeline_schedule',{'status':'Neighboring','percentage':2.0 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Neighboring


    ## neighbors
    if 'NB' in adata.uns['params']:
        adata = NB(adata)

    emit('get_pipeline_schedule',{'status':'Visualization','percentage':3.0 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Visualization

    ## 降维
    if 'DR' in adata.uns['params']:
        adata = DR(adata)
    else:
        adata.obsm['embedding'] = np.zeros((len(adata.obs),2))

    emit('get_pipeline_schedule',{'status':'Clustering','percentage':4.0 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Clustering

    ## 聚类
    if 'CL' in adata.uns['params']:
        adata = CL(adata)

    ## 保存聚类结果到localPoint
    ## 保存“供local pipeline进行划分时加载的缓存”
    localPointAdata.obs['label'] = adata.obs.label
    saveCache(localPointAdata, JobId, ViewId, 'localPoint')# uns为空，label只有


    # ## 过滤单样本的聚类，防止影响后续
    # adata = clearSimpleSizeCluster(adata)

    ## 轨迹推断
    if 'TI' in adata.uns['params']:
        adata = TI(adata)



    emit('get_pipeline_schedule',{'status':'Calculating DEGs','percentage':5.0 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Calculating Marker Genes

    ## 计算标志基因
    if 'MK' in adata.uns['params'] and len(adata.obs['label'].cat.categories) > 1:
        adata = MK(adata)
        ## 另外为raw_marker
        adata.uns['init_raw_marker'] = adata.uns['raw_marker']

    emit('get_pipeline_schedule',{'status':'Cell Chat','percentage':5.5/8},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Assembling Results

    ## 细胞通讯
    if 'CC' in adata.uns['params']:
        adata = CC(adata)

    emit('get_pipeline_schedule',{'status':'Assembling Results','percentage':6.0 / 7},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Assembling Results

    ## 保存 “最终缓存”
    saveCache(adata,adata.uns['JobId'],adata.uns['ViewId'],'Last')

    ## 构建保存metaData
    metaData = generateMetaDataFromAdata(adata)
    saveViewMetaData(JobId, ViewId, metaData)

    ## 构建保存Tree
    saveToTree(adata)


    return adata


def mergePipeline(params):
    adata = None

    return adata


'''
分析步骤
'''

## 采样 sampling
def SP(adata):
    '''
    adata: Anndata
    '''
    if 'Random' in adata.uns['params']['SP']:
        try:
            SPParams = {}
            if 'sampling_num' in adata.uns['params']['SP']['Random']:
                SPParams['sampling_num'] = adata.uns['params']['SP']['Random']['sampling_num']
            if 'sampling_radio' in adata.uns['params']['SP']['Random']:
                SPParams['sampling_radio'] = adata.uns['params']['SP']['Random']['sampling_radio']
            adata = sampleAdataRandom(adata,**SPParams)
            return adata
        except Exception as e:
            raise pipelineException(
                location='Sampling',
                advice='The current dataset cannot perform "Random Sampling", it may be because the data set contains too few cells.',
                message=str(e))
    elif 'Stratified' in adata.uns['params']['SP']:
        try:
            SPParams = {}
            if 'sampling_num' in adata.uns['params']['SP']['Stratified']:
                SPParams['sampling_num'] = adata.uns['params']['SP']['Stratified']['sampling_num']
            if 'sampling_radio' in adata.uns['params']['SP']['Stratified']:
                SPParams['sampling_radio'] = adata.uns['params']['SP']['Stratified']['sampling_radio']
            adata = sampleAdataStratified(adata,**SPParams)
            return adata
        except Exception as e:
            raise pipelineException(
                location='Sampling',
                advice='The current dataset cannot perform "Stratified Sampling", it may be because the data set contains too few cells.',
                message=str(e))

        
    

## 质量控制  quality control
def QC(adata):
    '''
    adata: Anndata
    '''
    ## 过滤掉异常细胞和异常基因，这里异常被定义为表达量低
    if 'filterCells' in adata.uns['params']['QC']:
        try:
            if 'min_genes' in adata.uns['params']['QC']['filterCells']:
                sc.pp.filter_cells(adata, min_genes = adata.uns['params']['QC']['filterCells']['min_genes'])
            if 'max_genes' in adata.uns['params']['QC']['filterCells']:
                sc.pp.filter_cells(adata, max_genes = adata.uns['params']['QC']['filterCells']['max_genes'])
        except Exception as e:
            raise pipelineException(
                location='Quality Control - Filter Outlying Cells',
                advice='The current dataset cannot perform "Filter Outlying Cells". Please either adjust the parameters for the "Filter Outlying Cells" or skip this step.',
                message=str(e))
    if 'filterGenes' in adata.uns['params']['QC']:
        try:
            if 'min_cells' in adata.uns['params']['QC']['filterGenes']:
                sc.pp.filter_genes(adata, min_cells = adata.uns['params']['QC']['filterGenes']['min_cells'])
            if 'max_cells' in adata.uns['params']['QC']['filterGenes']:
                sc.pp.filter_genes(adata, max_cells = adata.uns['params']['QC']['filterGenes']['max_cells'])
        except Exception as e:
            raise pipelineException(
                location='Quality Control - Filter Outlying Genes',
                advice='The current dataset cannot perform "Filter Outlying Genes". Please either adjust the parameters for the "Filter Outlying Genes" or skip this step.',
                message=str(e))
    ## 过滤掉高线粒体基因以及相关细胞
    if 'qcMetrics' in adata.uns['params']['QC']:
        try:
            organism = adata.uns['params']['QC']['qcMetrics']['type']
            if organism == 'Mouse':
                adata.var['mt'] = adata.var_names.str.startswith('mt-') #线粒体稳定
            else: ## HUMAN
                adata.var['mt'] = adata.var_names.str.startswith('MT-') #线粒体稳定
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            # if 'geneCounts' in adata.uns['params']['QC']['qcMetrics']:
            #     adata = adata[adata.obs.n_genes_by_counts < adata.uns['params']['QC']['qcMetrics']['geneCounts'], :]
            if 'pctCounts' in adata.uns['params']['QC']['qcMetrics']:
                adata = adata[adata.obs.pct_counts_mt < adata.uns['params']['QC']['qcMetrics']['pctCounts'], :]
            
        except Exception as e:
            raise pipelineException(
                location='Quality Control - MT-Gene-Control',
                advice='The current dataset cannot perform "MT-Gene-Control". Please either adjust the parameters for the "MT-Gene-Control" or skip this step.',
                message=str(e))

    return adata

## 数据集融合 Data Integration
def DI(adata):
    params = adata.uns['params']
    ### 编号batch
    batch = []
    batch.append(params['dataset']['name'])
    for item in adata.uns['params']['DI']['Datasets']:
        batch.append(item['Dataset']['name'])
    ### 取数据，并且进行简单的质量控制
    adatas = []
    for item in adata.uns['params']['DI']['Datasets']:
        tempData,hasNan = readData(item['Dataset'],adata.uns['params']['JobId'])
        tempData.uns['params'] = {}
        tempData.uns['params']['QC'] = item['qcParams']
        adatas.append(QC(tempData))
    ### 求基因的交集
    var_names = adata.var_names
    for item in adatas:
        var_names = var_names.intersection(item.var_names)
    adata = adata[:,var_names]
    adatas = list(map(lambda e:e[:,var_names],adatas))

    ## 融合
    ### ingest方法
    if 'Ingest' in adata.uns['params']['DI']['Method']:
        ### 处理ref数据集
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        if adata.uns['params']['DI']['Method']['Ingest']['embeddingMethod'] == 'UMAP':
            sc.tl.umap(adata)
        ### 处理query数据集
        if adata.uns['params']['DI']['Method']['Ingest']['embeddingMethod'] == 'UMAP':
            for i in range(len(adatas)):
                sc.tl.ingest(adatas[i],adata,embedding_method='umap')
        elif adata.uns['params']['DI']['Method']['Ingest']['embeddingMethod'] == 'PCA':
            for i in range(len(adatas)):
                sc.tl.ingest(adatas[i],adata,embedding_method='pca')
        ### 合并
        adata_all = adata.concatenate(*adatas,batch_categories=batch)
        if adata.uns['params']['DI']['Method']['Ingest']['embeddingMethod'] == 'UMAP':
            adata_all.obsm['di_prj'] = adata_all.obsm['X_umap']
        elif adata.uns['params']['DI']['Method']['Ingest']['embeddingMethod'] == 'PCA':
            adata_all.obsm['di_prj'] = adata_all.obsm['X_pca']
        adata_all.var['vst.variable'] = pd.Series(True,adata_all.var.index)
        ### 纠正
        adata_all = DIcorrect(adata_all)

    ### Scanorama方法
    # elif 'Scanorama' in adata.uns['params']['DI']['Method']:
    #     corrected = scanorama.correct_scanpy([adata,*adatas],return_dimred=True)
    #     adata_all = corrected[0].concatenate(*corrected[1:],batch_categories=batch)

    ### Harmony方法
    elif 'Harmony' in adata.uns['params']['DI']['Method']:
        adata_all = adata.concatenate(*adatas,batch_categories=batch)
        sc.tl.pca(adata_all)
        sc.external.pp.harmony_integrate(adata_all,'batch')
        adata_all.obsm['di_prj'] = adata_all.obsm['X_pca_harmony']
        adata_all.var['vst.variable'] = pd.Series(True,adata_all.var.index)
        ### 纠正
        adata_all = DIcorrect(adata_all)

    adata_all.uns['params'] = params
    return adata_all

## 变换 transform
def TS1(adata): 
    '''
    adata: Anndata
    '''
    if 'normalize' in adata.uns['params']['TS']:
        try:
            sc.pp.normalize_total(adata, target_sum=1e4)
        except Exception as e:
            raise pipelineException(
                location='Transformation - Normalize',
                advice='It could be that the dataset turned empty after previous preprocessing steps.',
                message=str(e))


    if 'log1p' in adata.uns['params']['TS']:
        ##校验输入数据中是否存在负数
        # if hasattr(adata.X,'A'):
        #     matrix = adata.X.A
        # else:
        #     matrix = adata.X
        matrix = get_dense_adata_X(adata.X)
        if (matrix < 0).any():#有inf存在
            raise pipelineException(
                location='Transformation - Log',
                advice='Please reconfigure the data preprocessing steps before applying the log function.',
                message='The data provided to the log function contains negative numbers, which is not valid.')
        try:
            sc.pp.log1p(adata)
        except Exception as e:
            raise pipelineException(
                location='Transformation - Log',
                advice='This could be because the dataset is empty at this moment.',
                message=str(e))


    # if 'CC' in adata.uns['params']:
    #     saveCache(adata,adata.uns['JobId'],adata.uns['ViewId'],'CellChat')


    return adata

def TS2(adata):
    '''
    adata:Anndata
    '''

    ## 当local pipeline或者使用了scTransform时，跳过该步骤
    ## TODO 为什么local就不能使用以下的两种变换方式？
    if 'local' not in adata.uns['params']['type'] and 'SCTransform' not in adata.uns['params']['FS']: 
        if 'regressOut' in adata.uns['params']['TS']:  
            try:  
                adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('MT-')
                sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) #TODO 这里是否和qc_metrics有强相关性？
            except Exception as e:
                raise pipelineException(
                    location='Transformation - Regress',
                    advice='No Advice',
                    message=str(e))

        if 'scale' in adata.uns['params']['TS']:
            try:
                sc.pp.scale(adata, max_value=10)
            except Exception as e:
                raise pipelineException(
                    location='Transformation - Scale',
                    advice='No Advice',
                    message=str(e))
    return adata


## 基因选择 Gene selection
def FS(adata):
    '''
    adata: Anndata
    '''
    if 'highlyVariableGenes' in adata.uns['params']['FS']:
        try:
            HVParams = {}   
            if 'topGenes' in adata.uns['params']['FS']['highlyVariableGenes']:
                HVParams['n_top_genes'] = adata.uns['params']['FS']['highlyVariableGenes']['topGenes']
            else:
                if 'minMean' in adata.uns['params']['FS']['highlyVariableGenes']:
                    HVParams['min_mean'] = adata.uns['params']['FS']['highlyVariableGenes']['minMean']
                else:
                    HVParams['min_mean'] = -np.inf
                if 'maxMean' in adata.uns['params']['FS']['highlyVariableGenes']:
                    HVParams['max_mean'] = adata.uns['params']['FS']['highlyVariableGenes']['maxMean']
                else:
                    HVParams['max_mean'] = np.inf
                if 'minDisp' in adata.uns['params']['FS']['highlyVariableGenes']:
                    HVParams['min_disp'] = adata.uns['params']['FS']['highlyVariableGenes']['minDisp']
                else:
                    HVParams['min_disp'] = -np.inf
                if 'maxDisp' in adata.uns['params']['FS']['highlyVariableGenes']:
                    HVParams['max_disp'] = adata.uns['params']['FS']['highlyVariableGenes']['maxDisp']
                else:
                    HVParams['max_disp'] = np.inf
            
            sc.pp.highly_variable_genes(adata, **HVParams)
            adata = adata[:, adata.var.highly_variable]
        except Exception as e:
            raise pipelineException(
                location='Gene Selection - HighlyVariable',
                advice='The error may be due to: 1. Incorrect parameter settings; 2. The current dataset does not have enough genes; 3. The current dataset contains invalid data. 4.The number of cells in a certain cluster in the data is too small.',
                message=str(e))
    elif 'scry' in adata.uns['params']['FS']:
        try:
            nTopGenes = adata.uns['params']['FS']['scry']['topGenes']
            # if hasattr(adata.uns['count'].X,'A'):
            #     adata = adata[:,scry(adata.uns['count'].X.A, adata.var.index, adata.obs.index,nTopGenes)]
            # else:
            #     adata = adata[:,scry(adata.uns['count'].X, adata.var.index, adata.obs.index,nTopGenes)]
            adata = adata[:,scry(get_dense_adata_X(adata.uns['count'].X), adata.var.index, adata.obs.index,nTopGenes)]
        except Exception as e:
            raise pipelineException(
                location='Gene Selection - scry',
                advice='The error may be due to: 1. Incorrect parameter settings; 2. The current dataset does not have enough genes; 3. The current dataset contains invalid data. 4.The number of cells in a certain cluster in the data is too small.',
                message=str(e))
    elif 'SCTransform' in adata.uns['params']['FS']:
        try:
            nTopGenes = adata.uns['params']['FS']['SCTransform']['topGenes']
            # if hasattr(adata.uns['count'].X,'A'):
            #     X = adata.uns['count'].X.A
            # else:
            #     X = adata.uns['count'].X
            X = get_dense_adata_X(adata.uns['count'].X)
            ## 将ndarray的数据类型转换为float32
            X = X.astype(np.float32)
            result = SCTransform(X, adata.var.index, adata.obs.index,nTopGenes,adata.uns['params']['JobId'])
            

            ##基因名词合法性检测 _ -> - 问题
            now_genes = result.columns.tolist()
            raw_genes = adata.var.index.tolist()
            for now_gene in now_genes:
                if now_gene not in raw_genes:
                    is_drop = True
                    for raw_gene in raw_genes:
                        if now_gene.replace("_","-") == raw_gene.replace("_","-"):#检测到匹配
                            is_drop = False
                            result.rename(columns={now_gene:raw_gene},inplace=True)
                            break
                        else:
                            pass
                    if is_drop: #丢弃基因
                        result.drop(now_gene,axis=1,inplace=True) #检查，是否列为基因


            ##过滤
            genes = result.columns.tolist()
            adata = adata[:,genes]
            cells = result.index.tolist()
            adata = adata[cells,:]
            # if hasattr(adata.uns['count'].X,'A'):
            if issparse(adata.uns['count'].X):
                adata.X = csr_matrix(result.to_numpy())
            else:
                adata.X = result.to_numpy()
            
        except Exception as e:
            raise pipelineException(
                location='Gene Selection - SCTransform',
                advice='The error may be due to: 1. The passed adata is not count; 2. The current dataset does not have enough genes; 3. The current dataset contains invalid data. 4.The number of cells in a certain cluster in the data is too small.',
                message=str(e))
    elif 'AllGenes' in adata.uns['params']['FS']:## 全基因选择
        pass
    elif 'marker' in adata.uns['params']['FS']:## only local mode
        try:
            if 'local' in adata.uns['params']['type']:
                ## Parent
                ParentId = adata.uns['params']['ParentId']
                JobId = adata.uns['params']['JobId']
                parentAdata = readCache(JobId, ParentId, 'Last')
                
                chosenCells = adata.obs.index.tolist()
                chosenLabels = parentAdata[chosenCells].obs.label.cat.categories.tolist()
                parentLabels = parentAdata.obs.label.cat.categories.tolist()
                parentMarkerNames = parentAdata.uns['raw_marker']['names']
                marker_list = []
                for label in chosenLabels:
                    i = parentLabels.index(label)
                    for j in range(len(parentMarkerNames)):
                        marker_list.append(parentMarkerNames[j][i])
                ## important: marker_list去重
                marker_list = list(set(marker_list))
                ## TODO 由于现在marker包括了adata last外的基因，所以需要需要进一步排除这些last外的基因 
                adata = adata[:,marker_list]
            elif 'combinedMarker' in  adata.uns['params']['FS']:## only local mode
                if 'local' in adata.uns['params']['type']:
                    ## Parent
                    ParentId = adata.uns['params']['ParentId']
                    JobId = adata.uns['params']['JobId']
                    parentAdata = readCache(JobId, ParentId, 'Last')
                    
                    chosenCells = adata.obs.index.tolist()
                    chosenLabels = parentAdata[chosenCells].obs.label.cat.categories.tolist()
                    parentLabels = parentAdata.obs.label.cat.categories.tolist()
                    parentMarkerNames = parentAdata.uns['raw_marker']['names']
                    marker_list = []
                    for label in chosenLabels:
                        i = parentLabels.index(label)
                        for j in range(len(parentMarkerNames)):
                            marker_list.append(parentMarkerNames[j][i])
                    ## important: marker_list去重
                    marker_list = list(set(marker_list))
                    ## TODO 由于现在marker包括了adata last外的基因，所以需要需要进一步排除这些last外的基因 

        except Exception as e:
            raise pipelineException(
                location='Gene Selection - marker',
                advice='The error may be due to: 1. Marker genes not being calculated; 2. The current dataset containing invalid data.',
                message=str(e))



    return adata

## 邻近图 neighbors
def NB(adata):
    '''
    adata: Anndata
    '''
    try:
        NBParams = {}
        sc.tl.pca(adata, svd_solver='arpack')
        if 'nNeighbors' in adata.uns['params']['NB']:
            NBParams['n_neighbors'] = adata.uns['params']['NB']['nNeighbors']
        if 'nPcs' in adata.uns['params']['NB']:
            NBParams['n_pcs'] = adata.uns['params']['NB']['nPcs']
        sc.pp.neighbors(adata, **NBParams)
    except Exception as e:
        raise pipelineException(
            location='Neighboring',
            advice='The error may be due to: an insufficient number of cells or the presence of invalid data in the dataset.',
            message=str(e))

    return adata

## 降维 dimension reduction
def DR(adata):  
    
    '''
    adata: Anndata
    '''
    ## 计算并保存TD
    try:
        TD = calculateNormTD(adata.X)
    except Exception as e:
        raise pipelineException(
            location='Visualization - calculate Distance Matrix',
            advice=':The error may be due to: 1. The dataset is too large, and there is insufficient memory to compute the distance matrix; 2. The current dataset contains invalid data."',
            message=str(e))

    ## 执行降维合并算法
    if 'UMAP' in adata.uns['params']['DR']:
        try:
            UMAP_params = {}
            if 'minDist' in adata.uns['params']['DR']['UMAP']:
                UMAP_params['min_dist'] = adata.uns['params']['DR']['UMAP']['minDist']
            if 'n_neighbors' in adata.uns['params']['DR']['UMAP']:
                UMAP_params['n_neighbors'] = adata.uns['params']['DR']['UMAP']['n_neighbors']
            
            ##是否预降维
            if 'PreDR' in adata.uns['params']['DR']['UMAP']:
                PreDR = adata.uns['params']['DR']['UMAP']['PreDR']
                if PreDR:
                    target_dimensions = 40
                    if adata.shape[1] > target_dimensions and adata.shape[0] > target_dimensions:
                        sc.pp.pca(adata,n_comps=target_dimensions) 
                    else:
                        adata.obsm['X_pca'] = adata.X
                    TD =  calculateNormTD(adata.obsm['X_pca'])
                    embedding = umap.UMAP(metric="precomputed",random_state= 0,**UMAP_params).fit_transform(TD)
                else:
                    embedding = umap.UMAP(metric="precomputed",random_state= 0,**UMAP_params).fit_transform(TD)
            else:
                embedding = umap.UMAP(metric="precomputed",random_state= 0,**UMAP_params).fit_transform(TD)
            adata.obsm['embedding'] = embedding
        except Exception as e:
            raise pipelineException(
                location='Visualization - UMAP',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data; 3. The parameter settings are unreasonable.',
                message=str(e))
    elif 'T-SNE' in adata.uns['params']['DR']:
        try:
            tSNE_params = {}
            if 'perplexity' in adata.uns['params']['DR']['T-SNE']:
                tSNE_params['perplexity'] = adata.uns['params']['DR']['T-SNE']['perplexity']
            ##是否预降维
            if 'PreDR' in adata.uns['params']['DR']['T-SNE']:
                PreDR = adata.uns['params']['DR']['T-SNE']['PreDR']
                if PreDR:
                    target_dimensions = 40
                    if adata.shape[1] > target_dimensions and adata.shape[0] > target_dimensions:
                        sc.pp.pca(adata,n_comps=target_dimensions) 
                    else:
                        adata.obsm['X_pca'] = adata.X
                    TD =  calculateNormTD(adata.obsm['X_pca'])
                    embedding = openTSNE.TSNE(metric="precomputed",random_state= 0,initialization='random',**tSNE_params).fit(TD)
                else:
                    embedding = openTSNE.TSNE(metric="precomputed",random_state= 0,initialization='random',**tSNE_params).fit(TD)
            else:
                embedding = openTSNE.TSNE(metric="precomputed",random_state= 0,initialization='random',**tSNE_params).fit(TD)
            embedding = np.array(embedding)
            adata.obsm['embedding'] = embedding
        except Exception as e:
            raise pipelineException(
                location='Visualization - t-SNE',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data; 3. The parameter settings are unreasonable.',
                message=str(e))
    elif 'PCA' in adata.uns['params']['DR']:
        try:
            sc.tl.pca(adata,n_comps=2,svd_solver='arpack')
            adata.obsm['embedding'] = adata.obsm['X_pca']
        except Exception as e:
            raise pipelineException(
                location='Visualization - PCA',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
                message=str(e))

    ## 保存TD
    adata.uns['TD'] = TD
    adata.uns['raw_TD'] = TD.copy()

    return adata

## 聚类 cluster
def CL(adata):
    '''
    adata: Anndata
    '''
    if 'leiden' in adata.uns['params']['CL']:
        try: #预降维
            target_dimensions = 40
            if adata.shape[1] > target_dimensions and adata.shape[0] > target_dimensions:
                sc.pp.pca(adata,n_comps=target_dimensions) 
            else:
                adata.obsm['X_pca'] = adata.X
                       
        except Exception as e:
            raise pipelineException(
                location='Clustering - leiden - Pre-dimensionality reduction',
                advice=':The error may be due to: 1. The dataset currently has too few genes or cells; 2. The current dataset contains invalid data;',
                message=str(e))
        try: ##neighbor
            neigh_params = {}
            if 'n_neighbors' in adata.uns['params']['CL']['leiden']:
                neigh_params['n_neighbors'] = adata.uns['params']['CL']['leiden']['n_neighbors']
                neigh_params['use_rep'] = 'X_pca'
            sc.pp.neighbors(adata,**neigh_params)
        except Exception as e:
            raise pipelineException(
                location='Clustering - leiden - neighboring',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
                message=str(e))
        try:
            leiden_params = {}
            if 'resolution' in adata.uns['params']['CL']['leiden']:
                leiden_params['resolution'] = adata.uns['params']['CL']['leiden']['resolution']
            sc.tl.leiden(adata, **leiden_params)
        except Exception as e:
            raise pipelineException(
                location='Clustering - leiden',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
                message=str(e))
        adata.obs['label'] = adata.obs['leiden']
    # elif 'louvain' in adata.uns['params']['CL']:
    #     try:
    #         neigh_params = {}
    #         if 'n_neighbors' in adata.uns['params']['CL']['louvain']:
    #             neigh_params['n_neighbors'] = adata.uns['params']['CL']['louvain']['n_neighbors']
    #             neigh_params['n_pcs'] = min(40,adata.shape[1])
    #         sc.pp.neighbors(adata, **neigh_params)
    #     except Exception as e:
    #         raise pipelineException(
    #             location='Clustering - louvain - neighboring',
    #             advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
    #             message=str(e))
    #     try:
    #         louvain_params = {}
    #         if 'resolution' in adata.uns['params']['CL']['louvain']:
    #             louvain_params['resolution'] = adata.uns['params']['CL']['louvain']['resolution']
    #         sc.tl.louvain(adata, **louvain_params)
    #     except Exception as e:
    #         raise pipelineException(
    #             location='Clustering - louvain',
    #             advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
    #             message=str(e))
    #     adata.obs['label'] = adata.obs['louvain']
    elif 'kmeans' in adata.uns['params']['CL']:
        try:
            ## 估计聚类数
            if 'auto_number' in adata.uns['params']['CL']['kmeans']:
                if adata.uns['params']['CL']['kmeans']['auto_number']:
                    X = get_dense_adata_X(adata.X)
                    ## 将ndarray的数据类型转换为float32
                    X = X.astype(np.float32)
                    recom_cluster_num = multik(X, adata.var.index, adata.obs.index,adata.uns['params']['JobId'])
        except Exception as e:
            raise pipelineException(
                location='Clustering - k-means - Auto Number',
                advice=':The error may be due to: 1. The dataset currently has too few genes or cells; 2. The current dataset contains invalid data;',
                message=str(e))

        try: # 预降维
            target_dimensions = 40
            if adata.shape[1] > target_dimensions and adata.shape[0] > target_dimensions:
                sc.pp.pca(adata,n_comps=target_dimensions) 
            else:
                adata.obsm['X_pca'] = adata.X          
        except Exception as e:
            raise pipelineException(
                location='Clustering - k-means - Pre-dimensionality reduction',
                advice=':The error may be due to: 1. The dataset currently has too few genes or cells; 2. The current dataset contains invalid data;',
                message=str(e))
        try: # 聚类
            kmeans_params = {}
            if 'n_clusters' in adata.uns['params']['CL']['kmeans']:
                kmeans_params['n_clusters'] = adata.uns['params']['CL']['kmeans']['n_clusters']
            if 'auto_number' in adata.uns['params']['CL']['kmeans']:
                if adata.uns['params']['CL']['kmeans']['auto_number']:
                    kmeans_params['n_clusters'] = int(recom_cluster_num)
            kmeans_params['random_state'] = 0           
            result = KMeans(**kmeans_params).fit(adata.obsm['X_pca']).labels_
            encoded_result = []
            for x in result:
                encoded_result.append(f'c_{x}')
            ## 转为categories的series
            categories = sorted(list(set(encoded_result)),key=lambda x:int(x[2:]))
            labels = pd.Series(pd.Categorical(encoded_result, categories=categories))
            labels.index = adata.obs.index
        except Exception as e:
            raise pipelineException(
                location='Clustering - k-means',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
                message=str(e))
        adata.obs['label'] = labels

    elif 'sc3s' in adata.uns['params']['CL']:
        try:
            ## 估计聚类数
            if 'auto_number' in adata.uns['params']['CL']['sc3s']:
                if adata.uns['params']['CL']['sc3s']['auto_number']:
                    X = get_dense_adata_X(adata.X)
                    ## 将ndarray的数据类型转换为float32
                    X = X.astype(np.float32)
                    recom_cluster_num = multik(X, adata.var.index, adata.obs.index,adata.uns['params']['JobId'])
        except Exception as e:
            raise pipelineException(
                location='Clustering - sc3s - Auto Number',
                advice=':The error may be due to: 1. The dataset currently has too few genes or cells; 2. The current dataset contains invalid data;',
                message=str(e))
            
        try: # 预降维
            target_dimensions = 40
            if adata.shape[1] > target_dimensions and adata.shape[0] > target_dimensions:
                sc.pp.pca(adata,n_comps=target_dimensions) 
            else:
                adata.obsm['X_pca'] = adata.X            
        except Exception as e:
            raise pipelineException(
                location='Clustering - sc3s - Pre-dimensionality reduction',
                advice=':The error may be due to: 1. The dataset currently has too few genes or cells; 2. The current dataset contains invalid data;',
                message=str(e))
        try: # 聚类
            sc3s_params = {}
            
            if 'n_clusters' in adata.uns['params']['CL']['sc3s']:
                final_cluster_num = adata.uns['params']['CL']['sc3s']['n_clusters']
            if 'auto_number' in adata.uns['params']['CL']['sc3s']:
                if adata.uns['params']['CL']['sc3s']['auto_number']:
                    final_cluster_num = int(recom_cluster_num)
            sc3s_params['n_clusters'] = final_cluster_num
            sc3s.tl.consensus(adata, **sc3s_params)
        except Exception as e:
            raise pipelineException(
                location='Clustering - sc3s',
                advice=':The error may be due to: 1. The dataset currently has too few cells; 2. The current dataset contains invalid data;',
                message=str(e))
        adata.obs['label'] = adata.obs['sc3s_' + str(final_cluster_num)]
        try: ##尾部处理
            adata.uns.pop('sc3s_trials')
        except:
            pass
    
    ## 合并单细胞聚类
    adata = mergeMiniSizeCluster(adata)
    
    ## refine the label name
    new_cate = []
    cluster_count = 0
    for item in adata.obs['label'].cat.categories:
        new_cate.append('c_' + str(cluster_count))
        cluster_count += 1
    adata.obs['label'] = adata.obs['label'].cat.rename_categories(new_cate)

    ## 聚类分层
    # sc.tl.dendrogram(adata,groupby='label',key_added='dendrogram')


    return adata

## 计算标志基因 marker gene
def MK(adata):
    '''
    adata: Anndata
    '''
    MarkerParams = {}
    if 'markerMethod' in adata.uns['params']['MK']:
        MarkerParams['method'] = adata.uns['params']['MK']['markerMethod']
    if 'nGenes' in adata.uns['params']['MK']:
        MarkerParams['n_genes'] = adata.uns['params']['MK']['nGenes']

    usedAdata = readCache(adata.uns['JobId'],adata.uns['ViewId'], 'Query')
    
    ## 过滤基因（确保usedAdata包含的基因都是adata中含有的）
    
    
    usedAdata.obs['label'] = adata.obs['label']

    if MarkerParams['method'] == 't-test': ## t-test
        try:
            sc.tl.rank_genes_groups(usedAdata,
                            groupby='label',
                            method = 't-test',
                            n_genes = MarkerParams['n_genes'],
                            tie_correct=False,
                            key_added='raw_marker')
        except Exception as e:
            raise pipelineException(
                location='Marker Identification - t-test',
                advice=':The error may be due to: 1. The dataset currently has too few genes; 2. The current dataset contains invalid data; 3. Clustering is not performed',
                message=str(e))
    elif MarkerParams['method'] == 't-test_overestim_var': ##  t-test(overestimate variance)
        try:
            sc.tl.rank_genes_groups(usedAdata,
                            groupby='label',
                            method = 't-test_overestim_var',
                            n_genes = MarkerParams['n_genes'],
                            tie_correct=False,
                            key_added='raw_marker')
        except Exception as e:
            raise pipelineException(
                location='Marker Identification - t-test(overestimate variance)',
                advice=':The error may be due to: 1. The dataset currently has too few genes; 2. The current dataset contains invalid data; 3. Clustering is not performed',
                message=str(e))
    elif MarkerParams['method'] == 'wilcoxon-test': ## wilcoxon-test
        try:
            sc.tl.rank_genes_groups(usedAdata,
                            groupby='label',
                            method = 'wilcoxon',
                            n_genes = MarkerParams['n_genes'],
                            tie_correct=False,
                            key_added='raw_marker')
        except Exception as e:
            raise pipelineException(
                location='Marker Identification - wilcoxon-test',
                advice=':The error may be due to: 1. The dataset currently has too few genes; 2. The current dataset contains invalid data; 3. Clustering is not performed',
                message=str(e))
    elif MarkerParams['method'] == 'wilcoxon-test(TLE)': ## wilcoxon-test(TLE)
        try:
            sc.tl.rank_genes_groups(usedAdata,
                            groupby='label',
                            method = 'wilcoxon',
                            n_genes = MarkerParams['n_genes'],
                            tie_correct=True,
                            key_added='raw_marker')
        except Exception as e:
            raise pipelineException(
                location='Marker Identification - wilcoxon-test(TLE)',
                advice=':The error may be due to: 1. The dataset currently has too few genes; 2. The current dataset contains invalid data; 3. Clustering is not performed',
                message=str(e))
    elif MarkerParams['method'] == 'logreg': ## logreg
        try:
            sc.tl.rank_genes_groups(usedAdata,
                            groupby='label',
                            method = 'logreg',
                            n_genes = MarkerParams['n_genes'],
                            key_added='raw_marker',
                            pts=True)
        except Exception as e:
            raise pipelineException(
                location='Marker Identification - logreg',
                advice=':The error may be due to: 1. The dataset currently has too few genes; 2. The current dataset contains invalid data; 3. Clustering is not performed',
                message=str(e))
    elif MarkerParams['method'] == 'COSG': ## COSG
        try:
            cosg.cosg(usedAdata,key_added='raw_marker',
                            mu=1,
                            n_genes_user=MarkerParams['n_genes'],
                            groupby='label')
        except Exception as e:
            raise pipelineException(
                location='Marker Identification - COSG',
                advice=':The error may be due to: 1. The dataset currently has too few genes; 2. The current dataset contains invalid data; 3. Clustering is not performed',
                message=str(e))

    ## 把raw marker附着
    adata.uns['raw_marker'] = usedAdata.uns['raw_marker']

    return adata

## 轨迹推断 Trajectory inference
def TI(adata):
    '''
    adata: Anndata
    '''
    if 'paga' in adata.uns['params']['TI']:
        sc.tl.paga(adata,groups='label')
        sc.pl.paga(adata,show=False,add_pos=True)
        sc.tl.draw_graph(adata,init_pos='paga')
    elif 'slingshot' in adata.uns['params']['TI']:
        rd = adata.obsm['embedding']
        cl = adata.obs['label'].to_list()
        result = slingshot(rd,cl)
        adata.uns['slingshot'] = result
    return adata

## 细胞通讯 Cell Chat
def CC(adata):
    '''
    adata:Anndata
    '''
    ## read data
    cc_data = readCache(adata.uns['JobId'],adata.uns['ViewId'],'Query')
    ## run
    if 'CellChat' in adata.uns['params']['CC']:
        try:
            # if hasattr(cc_data.X,'A'):
            #     cc_data_X = cc_data.X.A
            # else:
            #     cc_data_X = cc_data.X
            cc_data_X = get_dense_adata_X(cc_data.X)
            # result = CellChat(cc_data_X,cc_data.obs.index.tolist(),cc_data.var.index.tolist(),adata.obs['label'].tolist(),adata.uns['params']['CC']['CellChat']['organism'])
            result = CellChat(cc_data_X,cc_data.var.index,cc_data.obs.index,adata.obs['label'].tolist(),adata.uns['params']['CC']['CellChat']['organism'],adata.uns['params']['JobId'])
            adata.uns['CC'] = result
        except Exception as e:
            raise pipelineException(
                location='Cell Chat - Cell Chat',
                advice=':The error may be due to: 1. Clustering is not performed or Too few clusters; 2. Inproper dataset',
                message=str(e))
            
                    
    elif 'CellPhoneDB' in adata.uns['params']['CC']:
        try:
            # read lib
            CellChatDB = getCellChatDB(adata.uns['params']['CC']['CellPhoneDB']['curDB'][0],adata.uns['params']['CC']['CellPhoneDB']['curDB'][1])
            CellChatDB = CellChatDB.rename(columns={"ligand": "source", "receptor": "target"})
            
            # import dask
            # dask.config.set({"scheduler": "threads", "distributed.worker.daemon": False})
            # import sys
            # sys.modules["distributed"] = None

            ligrec_res = sq.gr.ligrec(
                adata,
                "label",
                interactions=CellChatDB,
                use_raw=False,
                copy=True
            )
            ligrec_res_mean = ligrec_res['means']

            cp_res = ligrec_res_mean.stack(level=['cluster_1', 'cluster_2']).reset_index()
            cp_res.columns = ['source', 'target', 'sender', 'receiver', 'value']
            cp_res = cp_res[cp_res['value'] != 0]


            ## 聚合
            source = []
            target = []
            ligand = []
            receptor = []
            prob = []
            score = []
            count_map = {}
            weight_map = {}
            for c_i in adata.obs['label']:
                count_map[c_i] = {}
                weight_map[c_i] = {}
                for c_j in adata.obs['label']:
                    count_map[c_i][c_j] = 0
                    weight_map[c_i][c_j] = 0
            for idx,row in cp_res.iterrows():
                source.append(row['sender'])
                target.append(row['receiver'])
                ligand.append(row['source'])
                receptor.append(row['target'])
                prob.append(row['value'])
                count_map[row['sender']][row['receiver']] += 1
                weight_map[row['sender']][row['receiver']] += row['value']
            adata.uns['CC'] = {
                'count':count_map,
                'weight':weight_map,
                'source':source,
                'target':target,
                'ligand':ligand,
                'receptor':receptor,
                'prob':prob
            }
        except Exception as e:
            raise pipelineException(
                location='Cell Chat - CellPhoneDB',
                advice=':The error may be due to: 1. Clustering is not performed or Too few clusters; 2. Inproper dataset',
                message=str(e))
    
    return adata


## 为adata生成新的metaData
def generateMetaDataFromAdata(adata):

    ## JobId
    JobId = adata.uns['JobId']
    ## ViewId
    ViewId = adata.uns['ViewId']
    ## ParentId
    ParentId = adata.uns['ParentId']


    ## generate
    metaData = {}
    group_ids = adata.obs.label.cat.categories
    ### color
    # colors = sns.color_palette("colorblind",n_colors=len(group_ids)).as_hex()
    colors = color_list[:len(group_ids)]
    metaData['group_color'] = {}
    for i,id in enumerate(group_ids):
        metaData['group_color'][id] = colors[i]
    ### name
    metaData['group_name'] = {}
    for i,id in enumerate(group_ids):
        metaData['group_name'][id] = id
    ### raw_embedding_range
    minRawEmbedding = np.min(adata.obsm['embedding'],axis=0)
    maxRawEmbedding = np.max(adata.obsm['embedding'],axis=0)
    metaData['raw_embedding_range'] = {
        'x':[float(minRawEmbedding[0]),float(maxRawEmbedding[0])],
        'y':[float(minRawEmbedding[1]),float(maxRawEmbedding[1])],
    }
    ### history_group_num
    metaData['history_group_num'] = len(group_ids)
    ### restore init
    metaData['init_state'] = {}
    metaData['init_state']['group_name'] = metaData['group_name']
    metaData['init_state']['group_color'] = metaData['group_color']
    metaData['init_state']['raw_embedding_range'] = metaData['raw_embedding_range']
    metaData['init_state']['history_group_num'] = metaData['history_group_num']
    if 'embedding' in adata.obsm:
        metaData['init_state']['raw_embedding'] = adata.obsm['embedding'].tolist()
    else:
        metaData['init_state']['raw_embedding'] = np.zeros((len(adata.obs),2)).tolist()
    metaData['init_state']['raw_labels'] = adata.obs['label'].tolist()



    return metaData


## 删除adata中单样本聚类
def clearSimpleSizeCluster(adata):
    labelIndex = {}
    for key,value in adata.obs.label.items():
        labelIndex[key] = (adata.obs['label'].value_counts() > 1)[value]
    labelIndex = pd.Series(labelIndex)
    return adata[labelIndex,:]

## 将adata中少样本的聚类合并到最近聚类中
def mergeMiniSizeCluster(adata):
    
    size_threshold = 3 #聚类合并的阈值 小于等于

    if 'label' not in adata.obs: ##如果没有标签，那么跳过
        return adata
    
    mini_clusters = (adata.obs['label'].value_counts() <= size_threshold).index[adata.obs['label'].value_counts() <= size_threshold].tolist()
    large_clusters = (adata.obs['label'].value_counts() > size_threshold).index[adata.obs['label'].value_counts() > size_threshold].tolist()
    
    if len(mini_clusters) == 0 or len(large_clusters) == 0:
        return adata
    
    min_cluster_adata = adata[adata.obs['label'].isin(mini_clusters)].copy()
    large_cluster_adata = adata[adata.obs['label'].isin(large_clusters)].copy()

    neigh = NearestNeighbors()
    neigh.fit(getXfromAdata(large_cluster_adata))
    nearest_indices = neigh.kneighbors(getXfromAdata(min_cluster_adata),n_neighbors=1,return_distance=False)

    new_min_cluster_labels = []
    for indice in range(len(nearest_indices)):
        new_min_cluster_labels.append(large_cluster_adata.obs['label'][nearest_indices[indice][0]])
    min_cluster_adata.obs['label'] = new_min_cluster_labels
    
    adata.obs['label'].update(min_cluster_adata.obs['label'])
    adata.obs['label'] = adata.obs['label'].cat.remove_unused_categories() ##去除不存在的类别
    
    removeUnusedLabelFromMetaData(adata)
    
    return adata