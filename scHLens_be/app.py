from cProfile import label
import enum
from unittest import result
from flask import Flask, request, jsonify,send_file,make_response,send_from_directory
import json
from flask.json import JSONEncoder
from zmq import DRAFT_API
from lib.pipeline import *
from lib.utils import *
from lib.vars import color_list
import scanpy as sc
import numpy as np
import pandas as pd
import umap
import time
import os
import config
import shutil
import traceback
import datetime
import zipfile
from flask_socketio import SocketIO, send,emit,join_room,leave_room,disconnect
from threading import Thread
from flask_cors import CORS
import locale
import datetime
from sklearn.preprocessing import MinMaxScaler
from types import MappingProxyType
import re
import asyncio

# 设置语言环境为英文
os.environ['LANG'] = 'en_US.UTF-8'
os.environ['LC_ALL'] = 'en_US.UTF-8'
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


app = Flask(__name__)
CORS(app)

socketio = SocketIO(app,cors_allowed_origins="*")


# 打开debug模式
app.debug = True
app.config.update(DEBUG=True)
app.config.from_object(config)
app.config.from_pyfile('config.cnf', silent=True)
app.json_encoder = NumpyEncoder


@app.route('/scHLens/api/',methods=['POST','GET'])
def hello_world(): 
    return 'Hello World!'





'''
关键
'''
## 流水线计算（global | local）
@app.route('/scHLens/api/runPipeline', methods=['POST'])
def runPipeline():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']

    try:
        if 'global' in reqParams['type']: ## 全局
            adata = globalPipeline(reqParams)
        elif 'local' in reqParams['type']: ## 局部
            adata = localPipeline(reqParams)

        emit('get_pipeline_schedule',{'status':'Returning Results','percentage':1},to=f'pipeline_schedule_{JobId}',namespace='/') ## pipeline percentage Info: Returning Results


    except pipelineException as e:
        traceback.print_exc()
        location = e.location
        advice = e.advice
        message = e.message
        response = {
            'type':'pipelineException',
            'attach':{
                'location':location,
                'advice':advice,
                'message':message,
            }
        }
        return make_response(response,500)

    except Exception as e:
        traceback.print_exc()
        response = {
            'type':'pipelineException',
            'attach':{
                'location':'unKnown',
                'advice':'unKnown',
                'message':str(e)
            }
        }
        return make_response(response,500)
    else:
        return {
            'JobId':adata.uns['JobId'],
            'ViewId':adata.uns['ViewId'],
            'status':1, ##成功
            'message':'success'
        }


## 从视图的adata中获取视图数据给前端
@app.route('/scHLens/api/fetchViewData', methods=['POST'])
def fetchViewData():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']

    adata = readCache(JobId,ViewId,'Last')

    dict_result = getResponseFromAdata(adata)

    return dict_result

## 删除视图
@app.route('/scHLens/api/requestDeleteView', methods=['POST'])
def requestDeleteView():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    deleteViewId = reqParams['deleteViewId']

    '''
    从tree中删除对应节点
    '''
    root = getTree(JobId)

    # 从树中查找要删除的节点（及其亲节点）
    def findNodeFromTree(cur,ToDeleteId): #ToDeleteId不能是根节点
        
        for child in cur['children']:

            if child['ViewId'] == ToDeleteId:
                return child,cur
            else:
                subResult,subParent = findNodeFromTree(child,ToDeleteId)
                if(subResult is not None):
                    return subResult,subParent
        return None,cur

    deletedNode,parent = findNodeFromTree(root,deleteViewId)

    # 从view文件中删除对应节点
    def deleteViewDataViaTree(cur):
        curViewId = cur['ViewId']
        dataPath = f'job/{JobId}/view/{curViewId}'
        metaPath = f'job/{JobId}/meta/{curViewId}.json'
        if os.path.exists(dataPath):
            shutil.rmtree(dataPath)
        if os.path.exists(metaPath):
            os.remove(metaPath)
        for child in cur['children']:
            deleteViewDataViaTree(child)
    deleteViewDataViaTree(deletedNode)

    #移除节点
    parent['children'].remove(deletedNode)
    
    saveTree(JobId,root)

    return 'success'


## 合并视图
@app.route('/scHLens/api/mergeViews',methods=['POST'])
def mergeViews():

    ## 读取数据
    reqParams = json.loads(request.get_data())

    try:

        JobId = reqParams['JobId']
        globalViewId = reqParams['globalViewId']
        localViewIdList = reqParams['localViewIdList']
        mergeOptions = reqParams['mergeOptions']
        isLabel = mergeOptions['isLabel']
        isProjection = mergeOptions['isProjection']
        isMergeSmallLabels = mergeOptions['isMergeSmallLabels']
        isMergeSameNameLabels = mergeOptions['isMergeSameNameLabels']


        globalAdata = readCache(JobId, globalViewId, 'Last')
        globalMetaData = getViewMetaData(JobId,globalViewId)


        '''
        embedding merge
        '''
        if isProjection:
            globalTD = globalAdata.uns['TD'] #这有一个默认记录的过程，改了globalTD相当于改了globalAdata.uns['TD']，因为是引用
            ## 距离矩阵替换
            for localViewId in localViewIdList:
                localAdata = readCache(JobId, localViewId, 'Last')
                localTD = localAdata.uns['TD']

                globalIndex = globalAdata.obs.index
                localIndex = localAdata.obs.index

                indexArr = [np.where(globalIndex==index)[0][0] for index in localIndex] # 对应的索引

                for i in range(len(indexArr)):
                    globalTD[indexArr[i],indexArr] = localTD[i] * 0.05
                
            ## 稳定性降维
            embedding = None
            if 'T-SNE' in globalAdata.uns['params']['DR']:
                tSNE_params = {}
                if 'perplexity' in globalAdata.uns['params']['DR']['T-SNE']:
                    tSNE_params['perplexity'] = globalAdata.uns['params']['DR']['T-SNE']['perplexity']
                embedding = openTSNE.TSNE(metric="precomputed",random_state= 0,initialization='random',**tSNE_params).fit(globalTD)
                embedding = np.array(embedding)
            elif 'UMAP' in globalAdata.uns['params']['DR']:
                UMAP_params = {}
                if 'minDist' in globalAdata.uns['params']['DR']['UMAP']:
                    UMAP_params['min_dist'] = globalAdata.uns['params']['DR']['UMAP']['minDist']
                if 'n_neighbors' in globalAdata.uns['params']['DR']['UMAP']:
                    UMAP_params['n_neighbors'] = globalAdata.uns['params']['DR']['UMAP']['n_neighbors']
                embedding = umap.UMAP(metric="precomputed", random_state= 0,**UMAP_params).fit_transform(globalTD)
                embedding = np.array(embedding)
            elif 'PCA' in globalAdata.uns['params']['DR']:
                embedding =  globalAdata.obsm['embedding'] #保留原始的embedding
            globalAdata.obsm['embedding'] = embedding

            ## set Raw Embedding Range
            minRawEmbedding = np.min(globalAdata.obsm['embedding'],axis=0)
            maxRawEmbedding = np.max(globalAdata.obsm['embedding'],axis=0)
            globalMetaData['raw_embedding_range'] = {
                'x':[float(minRawEmbedding[0]),float(maxRawEmbedding[0])],
                'y':[float(minRawEmbedding[1]),float(maxRawEmbedding[1])],
            }
        

        '''
        label adjust and merge labels
        '''
        if isLabel: ## cluster标签更新


            def __inner_removeUnusedLabelFromMetaData(adata,metaData):### 为metaData中删除不存在的类别
                
                if 'label' not in adata.obs:
                    return
                existLabels = adata.obs['label'].cat.categories.tolist()
                metaIds =  list(metaData['group_name'].keys())
                ## remove used group name
                for label in metaIds:
                    if label not in existLabels:
                        del metaData['group_name'][label]
                ## remove used group color
                for label in metaIds:
                    if label not in existLabels:
                        del metaData['group_color'][label]
                


            for localViewId in localViewIdList: ##把局部标签更新到全局标签上
                '''
                更新步骤：
                    对于每个localAdata，做到：
                        a. 给每个local label分配一个新id，用于global中
                        b. 给global metaData中新增对应的新id的属性（group name、group color）
                        c. 然后修改globalLabel中对应local的id
                    
                    要修改的量：
                        a. global metadata
                        b. globalAdata label value
                        c. globalAdata.cat.categories
                '''
                
                localAdata = readCache(JobId, localViewId, 'Last') ## 注意，localAdata千万不能保存，因为有做修改
                localMetaData = getViewMetaData(JobId,localViewId) ## 注意，localViewMetaData千万不能保存，因为有做修改
                                            
                IdMap = {}
                for oldlocalId in localAdata.obs['label'].cat.categories.tolist():## 对局部数据的每种类型的id进行遍历
                    ## 为local label -> global label 生成新id
                    newLocalId = f"c_{globalMetaData['history_group_num']}"
                    globalMetaData['history_group_num'] += 1
                    IdMap[oldlocalId] = newLocalId
                    ## add new id to global MetaData
                    globalMetaData['group_name'][newLocalId] = localMetaData['group_name'][oldlocalId] ### group name
                    exist_colors = globalMetaData['group_color'].values() ### group color
                    if localMetaData['group_color'][oldlocalId] not in exist_colors:
                        choose_color = localMetaData['group_color'][oldlocalId]
                    else:
                        choose_color = 'black'
                        for candidate_color in color_list:
                            if candidate_color not in exist_colors:
                                choose_color = candidate_color
                                break
                    globalMetaData['group_color'][newLocalId] = choose_color    
                ## 修改global中local对应的label
                newLocalAdataLabels = localAdata.obs['label'].cat.rename_categories(IdMap)
                newLocalAdataLabelsCat = newLocalAdataLabels.cat.categories.tolist()
                globalLabelsCat = globalAdata.obs['label'].cat.categories.tolist()
                globalAdata.obs['label'] = globalAdata.obs['label'].cat.add_categories(newLocalAdataLabelsCat)
                newLocalAdataLabels = newLocalAdataLabels.cat.add_categories(globalLabelsCat)
                globalAdata.obs['label'].update(newLocalAdataLabels)
                globalAdata.obs['label'] = globalAdata.obs['label'].cat.remove_unused_categories()
                ## 从global Metadata删除已经不存在的label
                __inner_removeUnusedLabelFromMetaData(globalAdata,globalMetaData)
            
            
            ##合并同名项
            if isMergeSameNameLabels: 
                globalIds =  globalAdata.obs['label'].cat.categories.tolist()
                for i in range(len(globalIds)):
                    for j in range(i+1,len(globalIds)):
                        id1 = globalIds[len(globalIds) - i -1]
                        id2 = globalIds[len(globalIds) - j -1]
                        if(globalMetaData['group_name'][id1] == globalMetaData['group_name'][id2]):###如果同名，id1合并到id2上
                            globalAdata.obs['label'] = globalAdata.obs['label'].replace(id1,id2)
                            globalAdata.obs['label'] = globalAdata.obs['label'].cat.remove_unused_categories()
                            break
                ## 从global Metadata删除已经不存在的label
                __inner_removeUnusedLabelFromMetaData(globalAdata,globalMetaData)
            
            ##合并小聚类            
            if isMergeSmallLabels: 
                size_threshold = 1 #聚类合并的阈值 小于等于
                mini_clusters = (globalAdata.obs['label'].value_counts() <= size_threshold).index[globalAdata.obs['label'].value_counts() <= size_threshold].tolist()
                large_clusters = (globalAdata.obs['label'].value_counts() > size_threshold).index[globalAdata.obs['label'].value_counts() > size_threshold].tolist()
                
                if len(mini_clusters) != 0 and len(large_clusters) != 0:
                    min_cluster_adata = globalAdata[globalAdata.obs['label'].isin(mini_clusters)].copy()
                    large_cluster_adata = globalAdata[globalAdata.obs['label'].isin(large_clusters)].copy()

                    neigh = NearestNeighbors()
                    neigh.fit(getXfromAdata(large_cluster_adata))
                    nearest_indices = neigh.kneighbors(getXfromAdata(min_cluster_adata),n_neighbors=1,return_distance=False)

                    new_min_cluster_labels = []
                    for indice in range(len(nearest_indices)):
                        new_min_cluster_labels.append(large_cluster_adata.obs['label'][nearest_indices[indice][0]])
                    min_cluster_adata.obs['label'] = new_min_cluster_labels
                    
                    globalAdata.obs['label'].update(min_cluster_adata.obs['label'])
                    globalAdata.obs['label'] = globalAdata.obs['label'].cat.remove_unused_categories() ##去除不存在的类别
                    
                    __inner_removeUnusedLabelFromMetaData(globalAdata,globalMetaData)
                    

            ## 重做marker
            if 'MK' in globalAdata.uns['params'] and len(globalAdata.obs['label'].cat.categories) > 1:
                globalAdata = MK(globalAdata)

        ## save changes
        saveViewMetaData(JobId,globalViewId,globalMetaData)
        saveCache(globalAdata,JobId,globalViewId,'Last')

    except Exception:
        traceback.print_exc()
        result = {
            'JobId':JobId,
            'ViewId':globalViewId,
            'status':0, ##报错
            'message':'error'

        }
    else:
        result = {
            'JobId':JobId,
            'ViewId':globalViewId,
            'status':1, ##成功
            'message':'success'
        }

    return result


## 恢复视图的原本标签（刚从pipeline出来）
@app.route('/scHLens/api/restoreViewLabels',methods=['POST'])
def restoreViewLabels():
    ## 读取参数
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']

    ## restore metaData
    metaData = getViewMetaData(JobId,ViewId)
    metaData['group_color'] = copy.deepcopy(metaData['init_state']['group_color'])
    metaData['group_name'] = copy.deepcopy(metaData['init_state']['group_name'])
    metaData['history_group_num'] = metaData['init_state']['history_group_num']

    ## restore adata
    adata = readCache(JobId, ViewId, 'Last')
    adata.obs['label'] = pd.Series(metaData['init_state']['raw_labels'],index=adata.obs.index.tolist()).astype('category')
    
    ## restore adata marker
    if 'init_raw_marker' in adata.uns:
        adata.uns['raw_marker'] = adata.uns['init_raw_marker']
    

    ## save
    saveViewMetaData(JobId,ViewId,metaData)
    saveCache(adata,JobId,ViewId,'Last')

    return 'success'




## 恢复视图的原本投影（刚从pipeline出来）
@app.route('/scHLens/api/restoreViewProjections',methods=['POST'])
def restoreViewProjections():
    ## 读取参数
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']

    ## restore metaData
    metaData = getViewMetaData(JobId,ViewId)
    metaData['raw_embedding_range'] = copy.deepcopy(metaData['init_state']['raw_embedding_range'])

    ## restore adata
    adata = readCache(JobId, ViewId, 'Last')
    adata.obsm['embedding'] = np.array(metaData['init_state']['raw_embedding'])
    adata.uns['TD'] = adata.uns['raw_TD']

    ## save
    saveViewMetaData(JobId,ViewId,metaData)
    saveCache(adata,JobId,ViewId,'Last')

    return 'success'

'''
Job
'''

## 创建Job
@app.route('/scHLens/api/createNewJob',methods=['POST'])
def createJob():
    JobId = initJob()
    return JobId

## 读取Job
@app.route('/scHLens/api/loadExistJob',methods=['POST'])
def loadJob():
    # 检索JobId是否存在
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    if os.path.exists('./job/'+ JobId):##JobId存在
        ##读取Tree信息
        tree = getTree(JobId)

        if not tree:
            return json.dumps({})

        def attachDataToTree(root):
            curView = root['ViewId']
            curData = readCache(JobId, curView, 'Last')
            curResult = getResponseFromAdata(curData)
            root['data'] = curResult
            if len(root['children']) != 0:
                for child in root['children']:
                    attachDataToTree(child)

        attachDataToTree(tree)
    
        ## 解决numpy数据无法进行json编码的问题
        class NpEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                                    np.int16, np.int32, np.int64, np.uint8,
                                    np.uint16, np.uint32, np.uint64)):
                    return int(obj)
                elif isinstance(obj, (np.float_, np.float16, np.float32,
                                    np.float64)):
                    return float(obj)
                elif isinstance(obj, (np.ndarray,)):
                    return obj.tolist()
                elif isinstance(obj, (np.bool_,)):
                    return bool(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        return json.dumps(tree,cls=NpEncoder)
    else:##JobId不存在
        return 'unexist'

## 导出job
@app.route('/scHLens/api/exportJob',methods=['POST'])
def exportJob():#TODO 修改时注意修改对应的回调函数

    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']

    ##压缩文件
    if not os.path.exists('./job_attachment/' + JobId):
        os.makedirs('job_attachment/' + JobId)
    shutil.make_archive(base_name='job_attachment/' + JobId + '/' + JobId, 
                        format='zip',
                        root_dir='./job',
                        base_dir='./' + JobId)


    path = './job_attachment/' + JobId + '/' + JobId +  '.zip'

    response = send_file(path,as_attachment=True)

    return response


## 上传job
@app.route('/scHLens/api/uploadJob',methods=['POST'])
def uploadJob():#TODO 修改时注意修改对应的回调函数

    def get_folder_names_in_zip(zip_path):#获取压缩包中顶层的文件夹名称
        folder_names = set()
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            for file_info in zip_ref.infolist():
                # 获取文件路径，并拆分成目录和文件名
                folder_path = file_info.filename
                if folder_path.endswith('/'):  # 检查路径是否以斜杠结尾，表示是文件夹
                    # 获取文件夹的顶级名称
                    top_level_folder = folder_path.split('/')[0]
                    folder_names.add(top_level_folder)
        return list(folder_names)

    _file = request.files.get('file')

    ## 保存压缩包
    archiveFilename = f'uploadJob_{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}.zip'
    archivePath = './globalAttachment/' + archiveFilename
    if not os.path.exists('globalAttachment'):
        os.makedirs('globalAttachment')
    _file.save(archivePath)


    uploadJobId = get_folder_names_in_zip(archivePath)[0]

    ## 如果job已存在，那么删除该job相关文件（覆盖操作）
    if os.path.exists('./job/' + uploadJobId):
        shutil.rmtree('./job/' + uploadJobId)
    if os.path.exists('./job_attachment/' + uploadJobId):
        shutil.rmtree('./job_attachment/' + uploadJobId)

    ## 解压压缩包
    shutil.unpack_archive(archivePath,'./job/')
    ## 创建job attachment文件夹
    os.makedirs('job_attachment/' + uploadJobId)

    ##清理上传压缩包
    if os.path.exists(archivePath):
        os.remove(archivePath)

    return jsonify({
        'uploadJobId':uploadJobId
    })


## 删除Job
@app.route('/scHLens/api/deleteJob',methods=['POST'])
def deleteJob():
    # 检索JobId是否存在
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    
    jobPath = f'job/{JobId}'
    jobAttachmentPath = f'job_attachment/{JobId}'
    if os.path.exists(jobPath):
        shutil.rmtree(jobPath)
    if os.path.exists(jobAttachmentPath):
        shutil.rmtree(jobAttachmentPath)


    return 'ok'

## 获取样本Job信息
@app.route('/scHLens/api/fetchSampleJobInfo',methods=['POST'])
def fetchSampleJobInfo():

    SampleJobInfo = []

    if not os.path.exists('sample_job'):#判断sample_job是否存在
        os.makedirs('sample_job')

    ##导入sample数据集
    for sample in os.listdir('./sample_job'):
        if os.path.exists('./sample_job/' + sample + '/config.json'):
            with open('./sample_job/' + sample + '/config.json','r',encoding = 'utf-8') as f:
                sample_config = json.load(f)
                SampleJobInfo.append(sample_config)

    return jsonify(SampleJobInfo)


## 加载Sample Job
@app.route('/scHLens/api/loadSampleJob',methods=['POST'])
def loadSampleJob():
    # 检索JobId是否存在
    reqParams = json.loads(request.get_data())
    sampleId = reqParams['sampleId']
    
    if not os.path.exists('sample_job'):#job sample
        os.makedirs('sample_job')
    
    ##根据现有sample创建新Job
    ### 获取新的JobId
    JobId = getNewJobId()
    
    ### 查找sample原本的路径
    samplePath = None
    for sample in os.listdir('./sample_job'):
        if os.path.exists('./sample_job/' + sample + '/config.json'):
            with open('./sample_job/' + sample + '/config.json','r',encoding = 'utf-8') as f:
                sample_config = json.load(f)
                if sample_config['id'] == sampleId:
                    samplePath = './sample_job/' + sample
                    break

    if samplePath is None:
        response = {
            'type':'cannotFindSampleException',
            'attach':{
                'location':'Load Sample',
                'advice':'Please Refresh the page, and then try again',
                'message':'cannot Find the Sample'
            }
        }
        return make_response(response,500)


    
    ### 创建job的文件夹
    JobPath = f'job/{JobId}'
    os.makedirs(JobPath)
    JobAttachmentPath = f'job_attachment/{JobId}'
    os.makedirs(JobAttachmentPath)
    ### 拷贝sample文件到job中
    for item in os.listdir(samplePath):
        source_item = f'{samplePath}/{item}'
        des_item = f'{JobPath}/{item}'
        if os.path.isfile(source_item): #文件
            shutil.copy2(source_item, des_item)
        elif os.path.isdir(source_item): #文件夹
            shutil.copytree(source_item, des_item)
            
            
    ### 把新job的中的adata cache的job id修改，要修改的有adata.uns['JobId']和adata.uns['params']['JobId']，adata.uns['count]中相应的JobId
    for view in os.listdir(f'{JobPath}/view'):
        viewPath = f'{JobPath}/view/{view}'
        for item in os.listdir(f'{viewPath}/cache'):
            cachePath = f'{viewPath}/cache/{item}'
            adata = sc.read(cachePath)
            ##解决存取时，log1p导致的小bug（对于字典中值为None的键，会直接在保存是忽略该键值）
            if 'log1p' in adata.uns and 'base' not in adata.uns['log1p']:
                adata.uns['log1p']['base'] = None
            ##修改
            adata.uns['JobId'] = JobId
            adata.uns['params']['JobId'] = JobId
            if 'count' in adata.uns:
                adata.uns['count'].uns['JobId'] = JobId
                adata.uns['count'].uns['params']['JobId'] = JobId

            adata.write(cachePath)

    
    return jsonify({
        'newJobId':JobId
    })

    

'''
数据集
'''

## 获取数据集信息
@app.route('/scHLens/api/fetchDatasets',methods=['POST'])
def fetchDatasets():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']

    datasetConfigs = []

    if not os.path.exists('sample'):#判断sample是否存在
        os.makedirs('sample')

    ##导入sample数据集
    for sample in os.listdir('./sample'):
        if os.path.exists('./sample/' + sample + '/config.json'):
            with open('./sample/' + sample + '/config.json','r',encoding = 'utf-8') as f:
                sample_config = json.load(f)
                datasetConfigs.append(sample_config)
    ##导入user上传的数据集
    for dataset in os.listdir('./job/' + JobId + '/dataset'):
        if os.path.exists('./job/' + JobId + '/dataset/' + dataset + '/config.json'):
            with open('./job/' + JobId + '/dataset/' + dataset + '/config.json','r',encoding = 'utf-8') as f:
                dataset_config = json.load(f)
                datasetConfigs.append(dataset_config)

    return jsonify(datasetConfigs)

## 上传数据集
@app.route('/scHLens/api/upload',methods=['POST']) #数据集文件上传函数，上传的每个文件都会调用一次该函数
def upload():
    upfile = request.files['file']
    dataset_name = request.form['name']
    dataset_type = request.form['type']
    JobId = request.form['JobId']
    filename = upfile.filename

    # 初始化文件夹
    if not os.path.exists('job/' + JobId + '/dataset/' + dataset_name):
        os.mkdir('job/' + JobId + '/dataset/' + dataset_name)

    # 装入数据集文件
    upfile.save('./job/' + JobId + '/dataset/' + dataset_name + '/' + filename)##保存本次上传的文件

    # 装入数据集配置文件
    if not os.path.exists('job/' + JobId + '/dataset/'  + dataset_name + '/config.json'):##如果之前没有config，那么装入config
        config = {}
        config['name'] = dataset_name
        config['type'] = dataset_type
        config['from'] = 'user'
        config_json = json.dumps(config)
        tempFile = open('job/' + JobId + '/dataset/' + dataset_name + '/config.json', 'w')
        tempFile.write(config_json)
        tempFile.close()
    
    return 'success'




## 检查上传的数据集是否有效，以及向config中补充dataset的shape
@app.route('/scHLens/api/check_dataset',methods=['POST']) #数据集文件上传函数，上传的每个文件都会调用一次该函数
def check_dataset():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    dataset_name = reqParams['dataset_name']
    
    dataset_config = getCustomDataSetConfig(JobId=JobId,dataset_name=dataset_name)
    try:
        adata,hasNan = readData(datasetConfig=dataset_config,JobId=JobId)
        if hasNan:
            emit('general_info',{
                'type':'warning',
                'title':'Invalid Values',
                'content':'This data set contains invalid values (such as null values).',
            },to=f'general_info_{JobId}',namespace='/')

    except Exception as e:
        ## remove invalid data set
        removeCustomDataSet(JobId=JobId,dataset_name=dataset_name)
        e.print_exc()
        ## 返回请求
        response = {
            'type':'invalidDatasetException',
            'attach':{
            }
        }
        return make_response(response,500)
    else:
        dataset_config['cell_num'] = adata.shape[0]
        dataset_config['gene_num'] = adata.shape[1]
        saveCustomDataSetConfig(JobId=JobId,dataset_name=dataset_name,dataset_config=dataset_config)
        
        return 'success'
        



## 上传gene sets
@app.route('/scHLens/api/uploadGeneSets',methods=['POST']) #gene set文件上传函数，上传的每个文件都会调用一次该函数
def uploadGeneSets():
    upfile = request.files['file']
    gene_set_name = request.form['name']
    gene_set_type = request.form['type']
    JobId = request.form['JobId']
    filename = f'set_{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}.json'

    ## 初始化type文件夹
    if not os.path.exists('job/' + JobId + '/geneset/' + gene_set_type):
        os.mkdir('job/' + JobId + '/geneset/' + gene_set_type)


    # 初始化具体gene sets文件夹
    if not os.path.exists('job/' + JobId + '/geneset/' + gene_set_type + '/' + gene_set_name):
        os.mkdir('job/' + JobId + '/geneset/' + gene_set_type + '/' + gene_set_name)


    # 装入数据集文件
    upfile.save('job/' + JobId + '/geneset/' + gene_set_type + '/' + gene_set_name + '/' + filename)##保存本次上传的文件
    
    # 在meta custom genesets中装入该数据集
    customGeneSetsConfig = getCustomGeneSetsConfig(JobId)
    if gene_set_type not in customGeneSetsConfig:
        customGeneSetsConfig[gene_set_type] = []
    customGeneSetsConfig[gene_set_type].append({
        'name':gene_set_name,
        'file':filename,
        'type':'custom'
    })
    setCustomGeneSetsConfig(JobId,customGeneSetsConfig)
    
    return 'success'


## 保存子数据集
@app.route('/scHLens/api/saveLocalDataset',methods=['POST'])
def saveLocalDataset(): #TODO 没有考虑数据融合的情况 #TODO 修改时注意修改对应的回调函数
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    chosenData = reqParams['chosenData']

    dataset,hasNan = readData(readCache(JobId, ViewId, 'Last').uns['params']['dataset'], JobId)
    
    ##过滤数据
    dataset = dataset[chosenData,:]

    ##导出为独立的数据集
    path = 'job/' + JobId + '/view/' + ViewId +  '/export/' + 'export' + '.h5ad' #TODO 这里的命名重复问题
    dataset.write(path)

    response = send_file(path,as_attachment=True,attachment_filename='export.h5ad')

    return response


'''
Group
'''

## 更改组名
@app.route('/scHLens/api/updateGroupName',methods=['POST'])
def updateGroupName():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    newGroupNames = reqParams['group_name']
    metaData = getViewMetaData(JobId, ViewId)
    metaData['group_name'] = newGroupNames
    saveViewMetaData(JobId, ViewId, metaData)

    return 'success'

## 更改颜色
@app.route('/scHLens/api/updateGroupColor',methods=['POST'])
def updateGroupColor():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    newGroupColors = reqParams['group_color']
    metaData = getViewMetaData(JobId, ViewId)
    metaData['group_color'] = newGroupColors
    saveViewMetaData(JobId, ViewId, metaData)

    return 'success'


## 合并重名标签
@app.route('/scHLens/api/mergeDuplicateLabels',methods=['POST'])
def mergeDuplicateLabels():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']

    adata = readCache(JobId=JobId,ViewId=ViewId,name='Last')
    metaData = getViewMetaData(JobId,ViewId)

    if 'label' not in adata.obs:
        return 'success'

    ##在obs.label中合并同名项
    Ids =  adata.obs['label'].cat.categories.tolist()
    for i in range(len(Ids)):
        for j in range(i+1,len(Ids)):
            id1 = Ids[len(Ids) - i -1]
            id2 = Ids[len(Ids) - j -1]
            if(metaData['group_name'][id1] == metaData['group_name'][id2]):###如果同名，id1合并到id2上
                adata.obs['label'] = adata.obs['label'].replace(id1,id2)
                adata.obs['label'] = adata.obs['label'].cat.remove_unused_categories()
                break
        
    ## 从global Metadata删除已经不存在的label
    existIds = adata.obs['label'].cat.categories.tolist()
    metaIds =  list(metaData['group_name'].keys())
    ### remove used group name
    for _id in metaIds:
        if _id not in existIds:
            del metaData['group_name'][_id]
    ### remove used group color
    for _id in metaIds:
        if _id not in existIds:
            del metaData['group_color'][_id]
    
    ## save changes
    saveViewMetaData(JobId,ViewId,metaData)
    saveCache(adata,JobId,ViewId,'Last')


    return 'success'

@app.route('/scHLens/api/exportGlobalMarkers',methods=['POST'])
def exportGlobalMarkers():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    GroupId = reqParams['GroupId']

    adata = readCache(JobId=JobId,ViewId=ViewId,name='Last')
    
    ## check marker(分为两种不同的情况讨论，global和local)
    if 'global' in adata.uns['params']['type']: # global模式
        if 'raw_marker' not in adata.uns:
            response = {
                'type':'MarkerUnexistException',
                'attach':{
                }
            }
            return make_response(response,500)
        raw_marker = adata.uns['raw_marker']
        cluster_id = GroupId
    else: # local模式
        rootViewId = getTree(JobId=JobId)['ViewId']
        globalAdata = readCache(JobId, rootViewId, 'Last')
        ## 构建新的cluster id
        cluster_id = '__TEMP__TARGET_CLUSTER'
        ## 替换标签
        globalLabels = globalAdata.obs['label'].astype('object')
        globalLabels.loc[adata[adata.obs['label'] == GroupId].obs.index]  = cluster_id
        globalAdata.obs['label'] = globalLabels.astype('category')
        ## 重做marker
        if 'MK' in globalAdata.uns['params'] and len(globalAdata.obs['label'].cat.categories) > 1 and len(globalAdata[globalAdata.obs['label']==cluster_id]) > 1: ## gloabl marker
            ### 过滤所有细胞数为1的类
            globalAdata = clearSimpleSizeCluster(globalAdata)
            globalAdata = MK(globalAdata)
            raw_marker = globalAdata.uns['raw_marker']
        else: ##不符合计算global的条件，使用局部marker
            if 'raw_marker' not in adata.uns:
                response = {
                    'type':'MarkerUnexistException',
                    'attach':{
                    }
                }
                return make_response(response,500)
            raw_marker = adata.uns['raw_marker']
            cluster_id = GroupId

    ## make rank
    gene_list = raw_marker['names'][cluster_id].tolist()

    ## save as excel
    df = pd.DataFrame(gene_list, columns=["Marker"])  
    path = 'job/' + JobId + '/view/' + ViewId +  '/export/' + 'exportGlobalMarkers' + '.xlsx' #TODO 这里的命名重复问题
    df.to_excel(path, index=False)

    response = send_file(path,as_attachment=True,attachment_filename='exportGlobalMarkers.xlsx')

    return response



@app.route('/scHLens/api/exportLocalMarkers',methods=['POST'])
def exportLocalMarkers():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    GroupId = reqParams['GroupId']

    adata = readCache(JobId=JobId,ViewId=ViewId,name='Last')
    
    ## 如果raw marker不存在，那么则人为做一个
    if 'raw_marker' not in adata.uns:
        response = {
            'type':'MarkerUnexistException',
            'attach':{
            }
        }
        return make_response(response,500)

    ## check marker
    raw_marker = adata.uns['raw_marker']
    cluster_id = GroupId

    ## get gene list
    gene_list = raw_marker['names'][cluster_id].tolist()

    ## save as excel
    df = pd.DataFrame(gene_list, columns=["Marker"])  
    path = 'job/' + JobId + '/view/' + ViewId +  '/export/' + 'exportLocalMarkers' + '.xlsx' #TODO 这里的命名重复问题
    df.to_excel(path, index=False)

    response = send_file(path,as_attachment=True,attachment_filename='exportLocalMarkers.xlsx')

    return response

'''
查询
'''

# 根据字符串查询符合匹配条件的基因数组
@app.route('/scHLens/api/queryCandidateGeneList', methods=['POST'])
def queryCandidateGeneList():
    reqParams = json.loads(request.get_data())
    ## read cache
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Query')
    geneList = adata.var.index[adata.var.index.str.contains(reqParams['geneMatch'], case=False)].tolist()
    return jsonify(geneList)

# 根据基因名查询该基因的表达值范围
@app.route('/scHLens/api/queryGeneValueRange', methods=['POST'])
def queryGeneValueRange():
    reqParams = json.loads(request.get_data())

    ## read cache
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Query')
    
    geneName = reqParams['geneName'] #要考虑可能是单个字符串基因，也可能是多个字符串组成的数组
    if not isinstance(geneName, str) and len(geneName) == 1:
        geneName = geneName[0]
    if isinstance(geneName,str):
        geneValueArr = None
        if hasattr(adata[:, [geneName]].X,'A'):
            geneValueArr = adata[:, [geneName]].X.A.flatten()
        else: 
            geneValueArr = adata[:, [geneName]].X.flatten()
    else:##包含多个基因
        sc.tl.score_genes(adata,gene_list=geneName)
        geneValueArr = adata.obs['score']
    
    return jsonify([float(geneValueArr.min()), float(geneValueArr.max())])

# 根据基因名查询每个细胞中的表达值数组
@app.route('/scHLens/api/queryGeneValueList', methods=['POST'])
def queryGeneValueList():
    reqParams = json.loads(request.get_data())

    ## read cache
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Query')

    keys = adata.obs_names.tolist()
    
    geneName = reqParams['geneName'] #要考虑可能是单个字符串基因，也可能是多个字符串组成的数组
    if not isinstance(geneName, str) and len(geneName) == 1:
        geneName = geneName[0]
    if isinstance(geneName,str):
        geneValueArr = None
        if hasattr(adata[:, [geneName]].X,'A'):
            geneValueArr = adata[:, [geneName]].X.A.flatten().tolist()
        else: 
            geneValueArr = adata[:, [geneName]].X.flatten().tolist()
    else:##包含多个基因
        sc.tl.score_genes(adata,gene_list=geneName)
        geneValueArr = adata.obs['score'].tolist()

    return jsonify(dict(zip(keys, geneValueArr)))

# 根据过滤器中的条件查询过滤细胞
@app.route('/scHLens/api/queryFilteredCellList', methods=['POST'])
def requestFilteredCellList():
    reqParams = json.loads(request.get_data())

    ## read cache
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Query')

    #adata = qualityControl(adata,reqParams)
    adata = adata[:, reqParams['geneName']]
    for index, item in enumerate(reqParams['geneRange']):
        if hasattr(adata.X,'A'):
            adata_X = adata.X.A
        else:
            adata_X = adata.X
        adata = adata[(adata_X[:, index] >= item[0]) & (adata_X[:, index] <= item[1])]
    return jsonify(adata.obs_names.tolist())


# 根据输入的多基因文本进行分割
@app.route('/scHLens/api/multiGenesSplitFromText', methods=['POST'])
def multiGenesSplitFromText():
    reqParams = json.loads(request.get_data())
    ## read cache
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Query')
    
    multiGeneText = reqParams['multiGeneText']
    split_genes = re.split(r'[\s;]+', multiGeneText.strip())
    # 过滤掉空字符串（如果有的话）
    split_genes = [item for item in split_genes if item]
    # 去重
    
    split_genes = list(set(split_genes))
    
    # 检测分割出的基因是否合法
    valid_genes = []
    invalid_genes = []
    ava_genes = adata.var.index.tolist()
    for gene in split_genes:
        if gene in ava_genes:
            valid_genes.append(gene)
        else:
            invalid_genes.append(gene)
    
    return jsonify({
        'split_num':len(split_genes),
        'valid_num':len(valid_genes),
        'invalid_num':len(invalid_genes),
        'valid_genes':valid_genes,
        'invalid_gens':invalid_genes
    })
        





'''
细胞类型推荐
'''
## 查询基因集信息
@app.route('/scHLens/api/queryGeneSets', methods=['POST'])
def queryGeneSets():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']

    gene_sets_info = getGeneSetsInfo(JobId)
    for org in gene_sets_info.keys():
        gene_sets_info[org] = list(map(lambda x:x['name'],gene_sets_info[org]))

    return jsonify(gene_sets_info)

## 查询gsea推荐结果(Prerank)
@app.route('/scHLens/api/queryGSEA', methods=['POST'])
def queryGSEA():
    reqParams = json.loads(request.get_data())

    ## read cache
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Last')

    ## read gene sets
    organism = reqParams['organism']
    gene_set_name = reqParams['gene_set_name']
    gene_sets = getGeneSet(JobId,organism,gene_set_name)

    ## 判断合法性
    if 'global' in adata.uns['params']['type'] and 'raw_marker' not in adata.uns: ## 如果全局没有marker可用
        return jsonify([])
    
    ## check marker(分为两种不同的情况讨论，global和local)
    if 'global' in adata.uns['params']['type']: # global模式
        raw_marker = adata.uns['raw_marker']
        cluster_id = reqParams['cluster_id']
    else: # local模式
        rootViewId = getTree(JobId=JobId)['ViewId']
        globalAdata = readCache(JobId, rootViewId, 'Last')
        ## 构建新的cluster id
        cluster_id = '__TEMP__TARGET_CLUSTER'
        ## 替换标签
        globalLabels = globalAdata.obs['label'].astype('object')
        globalLabels.loc[adata[adata.obs['label'] == reqParams['cluster_id']].obs.index]  = cluster_id
        globalAdata.obs['label'] = globalLabels.astype('category')
        ## 重做marker
        if 'MK' in globalAdata.uns['params'] and len(globalAdata.obs['label'].cat.categories) > 1 and len(globalAdata[globalAdata.obs['label']==cluster_id]) > 1:
            ### 过滤所有细胞数为1的类
            globalAdata = clearSimpleSizeCluster(globalAdata)
            globalAdata = MK(globalAdata)
            raw_marker = globalAdata.uns['raw_marker']
        else:
            raw_marker = adata.uns['raw_marker']
            cluster_id = reqParams['cluster_id'] 

    ## make rank
    gene_list = raw_marker['names'][cluster_id].tolist()
    value_list =raw_marker['scores'][cluster_id].tolist()
    logfoldchanges_threshold = reqParams['logfoldchanges_threshold']
    if 'logfoldchanges' in raw_marker and logfoldchanges_threshold != 'all': ## 按照logfoldchange进行过滤
        foldchange_list = raw_marker['logfoldchanges'][cluster_id].tolist()
        foldchange_list_remove_nan = np.where(np.isnan(foldchange_list), float('-inf'), foldchange_list)  ##防止nan值的影响
        new_gene_list = []
        new_value_list = []
        for i in range(0,len(gene_list)):
            if foldchange_list_remove_nan[i] > logfoldchanges_threshold:
                new_gene_list.append(gene_list[i])
                new_value_list.append(value_list[i])
        ## 防止logfoldchanges过滤得太狠，当总基因数小于10时，直接去前10个logfoldchanges值最高的基因及其对应值
        if len(new_gene_list) < 10:
            combined = zip(gene_list,value_list,foldchange_list_remove_nan)
            sorted_combined = sorted(combined, key=lambda x: x[2],reverse=True)
            sorted_gene_list,sorted_value_list,sorted_foldchange_list_remove_nan = zip(*sorted_combined) 
            new_gene_list = sorted_gene_list[:10]
            new_value_list = sorted_value_list[:10]
            ### 把得到的前10个基因，按照原本顺序在gene_list中的顺序排序
            combined2 = zip(new_gene_list,new_value_list)
            sorted_combined2 = sorted(combined2, key=lambda x: gene_list.index(x[0]))
            new_gene_list,new_value_list = zip(*sorted_combined2) 
        gene_list = new_gene_list
        value_list = new_value_list


    ## 排序gene_list（按照聚类在该基因上的平均表达值）,用以向用户推荐
    rootViewId = getTree(JobId=JobId)['ViewId'] # 获取当前job根节点id
    globalAdata = readCache(JobId, rootViewId, 'Query')
    def index_gene_value(gene):
        adata = globalAdata[globalAdata.obs['label'] == cluster_id,gene]
        if hasattr(adata.X,'A'):
            array = adata.X.A
        else:
            array = adata.X
        return np.mean(array)
    gene_list_SortedByExpression = sorted(gene_list,key=index_gene_value,reverse=True)

        
    if organism == 'Mouse':## 鼠基因转换
        gene_list,value_list = translateMouseGeneToHumanGene(gene_list,value_list)
    rnk = pd.Series(value_list,index=gene_list)

    ## run gsea preRank
    gsea_result = gp.prerank(rnk=rnk,gene_sets=gene_sets,min_size=1,max_size=len(gene_list)).res2d

    ## filter
    q_threshold = reqParams['q_threshold']
    p_threshold = reqParams['p_threshold']
    top = reqParams['top']
    gsea_result = gsea_result[gsea_result['ES'] > 0] ##按照ES > 0过滤
    if q_threshold != 'all':
        gsea_result = gsea_result[gsea_result['FDR q-val'] < q_threshold]##按照FDR q-val过滤
    if p_threshold != 'all':
        gsea_result = gsea_result[gsea_result['FWER p-val'] < p_threshold]##按照FWER p-val过滤

    ## sort
    gsea_result.sort_values(by=['FDR q-val','FWER p-val','NES'],ascending=[True,True,False],inplace=True)

    ##取前top个
    if top != 'all':
        gsea_result = gsea_result.head(top)
   

    return jsonify(gsea_result.to_dict(orient='records'))

## 查询enricher推荐结果
@app.route('/scHLens/api/queryEnricher', methods=['POST'])
def queryEnricher():
    reqParams = json.loads(request.get_data())

    ## read cache
    JobId = reqParams['JobId'] 
    ViewId = reqParams['ViewId']
    adata = readCache(JobId, ViewId, 'Last')

    ## read gene sets
    organism = reqParams['organism']
    gene_set_name = reqParams['gene_set_name']
    gene_sets = getGeneSet(JobId,organism,gene_set_name)


    ## 判断合法性
    if 'global' in adata.uns['params']['type'] and 'raw_marker' not in adata.uns: ## 如果全局没有marker可用
        return jsonify([])
    
    ## check marker(分为两种不同的情况讨论，global和local)
    if 'global' in adata.uns['params']['type']: # global模式
        raw_marker = adata.uns['raw_marker']
        cluster_id = reqParams['cluster_id']
    else: # local模式
        rootViewId = getTree(JobId=JobId)['ViewId']
        globalAdata = readCache(JobId, rootViewId, 'Last')
        ## 构建新的cluster id
        cluster_id = '__TEMP__TARGET_CLUSTER'
        ## 替换标签
        globalLabels = globalAdata.obs['label'].astype('object')
        globalLabels.loc[adata[adata.obs['label'] == reqParams['cluster_id']].obs.index]  = cluster_id
        globalAdata.obs['label'] = globalLabels.astype('category')
        ## 重做marker
        if 'MK' in globalAdata.uns['params'] and len(globalAdata.obs['label'].cat.categories) > 1 and len(globalAdata[globalAdata.obs['label']==cluster_id]) > 1:
            ### 过滤所有细胞数为1的类
            globalAdata = clearSimpleSizeCluster(globalAdata)
            globalAdata = MK(globalAdata)
            raw_marker = globalAdata.uns['raw_marker']
        else:
            raw_marker = adata.uns['raw_marker']
            cluster_id = reqParams['cluster_id']

    ## get gene_list（为了省事，这里连着value_list一起计算了）
    gene_list = raw_marker['names'][cluster_id].tolist()
    value_list = raw_marker['scores'][cluster_id].tolist()
    logfoldchanges_threshold = reqParams['logfoldchanges_threshold']
    if 'logfoldchanges' in raw_marker and logfoldchanges_threshold != 'all': ## 按照logfoldchange进行过滤
        foldchange_list = raw_marker['logfoldchanges'][cluster_id].tolist()
        foldchange_list_remove_nan = np.where(np.isnan(foldchange_list), float('-inf'), foldchange_list)  ##防止nan值的影响
        new_gene_list = []
        new_value_list = []
        for i in range(0,len(gene_list)):
            if foldchange_list_remove_nan[i] > logfoldchanges_threshold:
                new_gene_list.append(gene_list[i])
                new_value_list.append(value_list[i])
        ## 防止logfoldchanges过滤得太狠，当总基因数小于10时，直接去前10个logfoldchanges值最高的基因及其对应值
        if len(new_gene_list) < 20:
            combined = zip(gene_list,value_list,foldchange_list_remove_nan)
            sorted_combined = sorted(combined, key=lambda x: x[2],reverse=True)
            sorted_gene_list,sorted_value_list,sorted_foldchange_list_remove_nan = zip(*sorted_combined) 
            new_gene_list = sorted_gene_list[:20]
            new_value_list = sorted_value_list[:20]
            ### 把得到的前12个基因，按照原本顺序在gene_list中的顺序排序
            combined2 = zip(new_gene_list,new_value_list)
            sorted_combined2 = sorted(combined2, key=lambda x: gene_list.index(x[0]))
            new_gene_list,new_value_list = zip(*sorted_combined2) 
        gene_list = new_gene_list
        value_list = new_value_list


    ## 排序gene_list（按照聚类在该基因上的平均表达值）,用以向用户推荐
    rootViewId = getTree(JobId=JobId)['ViewId'] # 获取当前job根节点id
    globalAdata = readCache(JobId, rootViewId, 'Query')
    def index_gene_value(gene):
        adata = globalAdata[globalAdata.obs['label'] == cluster_id,gene]
        if hasattr(adata.X,'A'):
            array = adata.X.A
        else:
            array = adata.X
        return np.mean(array)
    gene_list_SortedByExpression = sorted(gene_list,key=index_gene_value,reverse=True)
        
    if organism == 'Mouse':## 鼠基因转换
        gene_list,value_list = translateMouseGeneToHumanGene(gene_list,value_list)




    ## run enrichr
    enrichr_result = gp.enrichr(gene_list,
                        gene_sets=gene_sets,
                        outdir=None).res2d

    ## filter
    p_threshold = reqParams['p_threshold']
    top = reqParams['top']

    if p_threshold != 'all':##按照Adjusted P-value过滤
        enrichr_result = enrichr_result[enrichr_result['Adjusted P-value'] < p_threshold]

    ## sort
    enrichr_result.sort_values(by=['Adjusted P-value','Combined Score'],ascending=[True,False],inplace=True)

    ##取前top个
    if top != 'all':
        enrichr_result = enrichr_result.head(top)


    return jsonify(enrichr_result.to_dict(orient='records'))




'''
其他
'''
# # 基因推荐（marker）
# @app.route('/scHLens/api/recommendGene',methods=['POST'])
# def recommendGene():
#     reqParams =  json.loads(request.get_data())

#     ## read cache
#     JobId = reqParams['JobId']
#     ViewId = reqParams['ViewId']
#     adata = readCache(JobId, ViewId, 'Last')
    
#     if reqParams['mode'] == 'HighlyVariable':
#         return jsonify(adata[:,adata.var.highly_variable].var.index.tolist())
#     elif reqParams['mode'] == 'Marker':
#         np.unique(adata.obs['label']).tolist()
#         return jsonify([])

## 用户打开的实例页面被关闭的事件触发函数
@app.route('/scHLens/api/InstanceClose', methods=['POST'])
def InstanceClose():
    reqParams = json.loads(request.get_data())
    
    return jsonify('close')

## 删除选择的节点
@app.route('/scHLens/api/updateDeleteCells', methods=['POST'])
def updateDeleteCells():
    reqParams = json.loads(request.get_data())
    JobId = reqParams['JobId']
    ViewId = reqParams['ViewId']
    chosenData = reqParams['chosenData']

    ## 更新deleteCells
    deleteCells = getDeleteCells(JobId)
    newDeleteCells = list(set(deleteCells) | set(chosenData))
    setDeleteCells(JobId,newDeleteCells)

    ## 视图更新Tree
    EmbeddingTree = getTree(JobId)
    def addEmbeddingToTree(cur):
        adata = readCache(JobId,cur['ViewId'],'Last')
        newResult = getResponseFromAdata(adata) ##由于deleteCells更新，所以从adata提取的数据也要取最新的
        cur['updateData'] = {
                'cellData':newResult['cellData'],
                'groups':newResult['groups'],
                'MK':newResult['MK'],
                'globalScores':newResult['globalScores'],
                'localScores':newResult['localScores'],
            }
        for child in cur['children']:
            addEmbeddingToTree(child)
    addEmbeddingToTree(EmbeddingTree)
    

    return jsonify(EmbeddingTree)


## 提交意见
@app.route('/scHLens/api/sendMessage', methods=['POST'])
def sendMessage():
    form = json.loads(request.get_data())
    first_name = form['first_name']
    last_name = form['last_name']
    email = form['email']
    message = form['message']
    
    if not os.path.exists('message'):
        os.makedirs('message')
    
    now_time = datetime.datetime.now()
    ## 四位随机ID
    characterList = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    ID = ''.join(random.choice(characterList) for _ in range(4))
    filename = now_time.strftime("%Y%m%d%H%M%S") + '-' + ID + '.txt'
    
    store_str = \
        'First name: ' + first_name + '\n' + \
        'Last name: ' + last_name + '\n' + \
        'Email: ' + email + '\n' + \
        'Time: ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + '\n' + \
        'Message: ' + message
    
    with open('message/' + filename,'w',encoding='utf-8') as f:
        f.write(store_str)
    return jsonify(1)
        


'''
socket 监听
'''
@socketio.on('connect')
def handle_connect(auth):
    print('socket connected!')

    if auth['type'] == 'Test':# 检测客户端是否启动
        emit('test',{'status':0}) #表示客户端已经成功启动
        # disconnect()
    elif auth['type'] == 'Main':
        JobId = auth['JobId']
        join_room(f'pipeline_schedule_{JobId}')
        join_room(f'general_info_{JobId}')

    return

@socketio.on('disconnect')
def handle_disconnect():
    print('socket disconnected!')
    return


@socketio.on('message')
def handle_message(message):
    return


'''
请求回调函数
'''
@app.route('/scHLens/api/clearCallback',methods=['POST'])
def clearCallback():
    reqParams =  json.loads(request.get_data())
    JobId = reqParams['JobId']
    callbackUrl = reqParams['url']
    params = reqParams['params']
    if callbackUrl == 'api/exportJob':# 导出Job
        exportPath = './job_attachment/' + JobId + '/' + JobId +  '.zip'
        if os.path.exists(exportPath):
            os.remove(exportPath)
    elif callbackUrl == 'api/saveLocalDataset': # 导出选中细胞作为数据集
        ViewId = params['ViewId']
        exportPath = 'job/' + JobId + '/view/' + ViewId +  '/export/' + 'export' + '.h5ad'       
        if os.path.exists(exportPath):
            os.remove(exportPath)
    elif callbackUrl == 'api/exportLocalMarkers': # 导出Local Markers
        ViewId = params['ViewId']
        exportPath = 'job/' + JobId + '/view/' + ViewId +  '/export/' + 'exportLocalMarkers' + '.xlsx'
        if os.path.exists(exportPath):
            os.remove(exportPath)
    elif callbackUrl == 'api/exportGlobalMarkers': # 导出Global Markers
        ViewId = params['ViewId']
        exportPath = 'job/' + JobId + '/view/' + ViewId +  '/export/' + 'exportGlobalMarkers' + '.xlsx'
        if os.path.exists(exportPath):
            os.remove(exportPath)
    return 'ok'

if __name__ == '__main__': ##!important vscode的debug不会走该路径执行该函数..
    ## 初始化
    initApp() 

    # app.run(port = 5003,debug=True)
    
    socketio.run(app,port=5003,debug=True,allow_unsafe_werkzeug=True)


