import os

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.sparse as sps
import scanpy as sc
# import scanorama
from anndata import AnnData

##校正合并后的矩阵
'''
adata_all:合并后的adata,带有.obsm['di_prj']、.obs['batch']以及.var['vst.variable']
'''
def DIcorrect(adata_all):

    datasets = []  
    datasets_dimred = []   # 存 harmony 的 embedding
    adatas = []
    batchs = adata_all.obs['batch'].unique()
    genes = adata_all.var.index

    for b in batchs:
        cells = adata_all.obs['batch'] == b   # 按样本选数据
        temp = adata_all[cells,adata_all.var['vst.variable']].copy()
        temp.X = csr_matrix(temp.X)
        
        adatas.append(temp)
        datasets.append(temp.X)
        datasets_dimred.append(temp.obsm['di_prj'])

    ALPHA = 0.10
    APPROX = True
    BATCH_SIZE = 5000
    DIMRED = 100
    HVG = None
    KNN = 20
    N_ITER = 500
    PERPLEXITY = 1200
    SIGMA = 15
    VERBOSE = 2
    ds_names = None

    # datasets_dimred = scanorama.assemble(datasets_dimred, 
    #             expr_datasets=datasets,   # Modified in place.
    #                                     verbose=VERBOSE, knn=KNN, sigma=SIGMA, approx=APPROX,
    #                                     alpha=ALPHA, ds_names=ds_names, batch_size=BATCH_SIZE,
    #                                     )

    new_adatas = []
    genes = adata_all[:, adata_all.var['vst.variable']].var_names.to_list()
    for i in range(len((adatas))):
        adata = AnnData(datasets[i])
        adata.obs = adatas[i].obs
        adata.obsm = adatas[i].obsm
        
        adata.var_names = genes

        adata.uns = adatas[i].uns
        new_adatas.append(adata)


    adata_corr = sc.concat(new_adatas)

    return adata_corr
