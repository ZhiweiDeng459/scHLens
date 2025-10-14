import axios from "axios";
import io from 'socket.io-client'

/**
 * 
 * 
 * 前端与后端的网络接口
 * 
 * 
 */



/**
 * 
 * 关键
 * 
 */

//流水线计算
export function runPipeline(params){
    return axios({
        method:"post",
        url:"api/runPipeline",
        data:params
    })
}

//获取单个视图的数据
export function fetchViewData(JobId,ViewId){
    return axios({
        method:"post",
        url:"api/fetchViewData",
        data:{
            JobId,
            ViewId
        }
    })

}

//删除视图
export function requestDeleteView(JobId,deleteViewId){
    return axios({
        method:"post",
        url:"api/requestDeleteView",
        data:{
            'JobId':JobId,
            'deleteViewId':deleteViewId,
        },
    })
}


//视图合并 值得注意的是，localPerIdList在后面的ID的有替换高优先级（即最后替换，可以覆盖前面的）
export function mergeViews(JobId,globalViewId,localViewIdList,mergeOptions){
    return axios({
        method:"post",
        url:"api/mergeViews",
        data:{
            'JobId':JobId,
            'globalViewId':globalViewId,
            'localViewIdList':localViewIdList,
            'mergeOptions':mergeOptions,
        },
    })
}


//恢复投影
export function restoreViewProjections(JobId,ViewId){
    return axios({
        method:"post",
        url:"api/restoreViewProjections",
        data:{
            'JobId':JobId,
            'ViewId':ViewId,
        },
    })
}

//恢复标签
export function restoreViewLabels(JobId,ViewId){
    return axios({
        method:"post",
        url:"api/restoreViewLabels",
        data:{
            'JobId':JobId,
            'ViewId':ViewId,
        },
    })
}


/**
 * 
 * Job
 * 
 */

//创建新的job
export function createNewJob(){
    return axios({
        method:"post",
        url:"api/createNewJob",
    })
}

//获取旧的job信息
export function loadExistJob(JobId){
    return axios({
        method:"post",
        url:"api/loadExistJob",
        data:{
            'JobId':JobId
        }
    })
}

//导出job
export function exportJob(JobId,handleProcess){
    return axios({
        method:"post",
        url:"api/exportJob",
        responseType: 'blob',//对于传输文件很重要，一定要添加，不然文件损坏
        data:{
            'JobId':JobId
        },
        onDownloadProgress:handleProcess
    })
}


//删除job
export function deleteJob(JobId){
    return axios({
        method:"post",
        url:"api/deleteJob",
        data:{
            'JobId':JobId
        }
    })
}



//上传job
export function uploadJob(formData){
    return axios({
        url:"api/uploadJob",
        method:"post",
        data:formData,
        'Content-type' : 'multipart/form-data',
    })
}


//清理数据的回调函数（通常在exportJob,saveLocalDataset之后调用）
export function clearCallback(JobId,url,params={}){
    return axios({
        method:"post",
        url:"api/clearCallback",
        data:{
            'JobId':JobId,
            'url':url,
            'params':params
        }
    })
}

//获取样本Job
export function fetchSampleJobInfo(){
    return axios({
        method:"post",
        url:"api/fetchSampleJobInfo",
        data:{
        },        
    })
}

export function loadSampleJob(sampleId){//加载SampleJob
    return axios({
        method:"post",
        url:"api/loadSampleJob",
        data:{
            'sampleId':sampleId
        },        
    })
}

/**
 * 
 * DataSet
 *  
 */


//获取数据集
export function fetchDatasets(JobId){
    return axios({
        method:"post",
        url:"api/fetchDatasets",
        data:{
            'JobId':JobId
        },        
    })
}


//保存子数据集
export function saveLocalDataset(JobId,ViewId,chosenData){
    return axios({
        method: "post",
        url: "api/saveLocalDataset",
        responseType: 'blob',
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'chosenData':chosenData, 
        },
    });
}



/**
 * 
 * Group
 * 
 */
//更改组名
export function updateGroupName(JobId,ViewId,group_name){
    return axios({
        method: "post",
        url: "api/updateGroupName",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'group_name':group_name, 
        },
    });   
}

//更改颜色
export function updateGroupColor(JobId,ViewId,group_color){
    return axios({
        method: "post",
        url: "api/updateGroupColor",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'group_color':group_color, 
        },
    });
}

//合并同名标签
export function mergeDuplicateLabels(JobId,ViewId){
    return axios({
        method: "post",
        url: "api/mergeDuplicateLabels",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
        },
    });
}


//导出当前组的global marker基因
export function exportGlobalMarkers(JobId,ViewId,GroupId){
    return axios({
        method: "post",
        url: "api/exportGlobalMarkers",
        responseType: 'blob',
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'GroupId':GroupId,
        },
    });

}

//导出当前组的local marker基因
export function exportLocalMarkers(JobId,ViewId,GroupId){
    return axios({
        method: "post",
        url: "api/exportLocalMarkers",
        responseType: 'blob',
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'GroupId':GroupId,
        },
    });

}



/**
 * 
 * 查询
 * 
 */

// 根据字符串查询符合匹配条件（包含）的基因数组
export function requestCandidateGeneList(JobId,ViewId,geneMatch) {
    return axios({
        method: "post",
        url: "api/queryCandidateGeneList",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'geneMatch':geneMatch, 
        },
    });
}

// 根据基因名查询该基因的表达值范围
export function requestGeneValueRange(JobId,ViewId,geneName) {
    return axios({
        method: "post",
        url: "api/queryGeneValueRange",
        data: { 
            'JobId':JobId,
            'ViewId':ViewId,
            'geneName':geneName,//可能是一个gene字符串，也可能是一个gene字符串数组
        },
    });
}

// 根据基因名查询每个细胞中的表达值数组
export function requestGeneValueList(JobId,ViewId,geneName) {
    return axios({
        method: "post",
        url: "api/queryGeneValueList",
        data: { 
            'JobId':JobId,
            'ViewId':ViewId,
            'geneName' : geneName //可能是一个gene字符串，也可能是一个gene字符串数组
        },
    });
}

// 根据基因名和表达值范围查询符合筛选条件的细胞数组
export function requestFilteredCellList(JobId,ViewId,geneName, geneRange) {
    return axios({
        method: "post",
        url: "api/queryFilteredCellList",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'geneName':geneName,
            'geneRange':geneRange, 
        },
    });
}

//根据输入的多基因文本进行分割
export function multiGenesSplitFromText(JobId,ViewId,multiGeneText) {
    return axios({
        method: "post",
        url: "api/multiGenesSplitFromText",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'multiGeneText':multiGeneText
        },
    });
}



/**
 * 
 * 细胞通讯
 * 
 */

export function queryCellChatDB(){
    return axios({
        method: "post",
        url: "api/queryCellChatDB",
        data: {
        },
    })
}

/**
 * 
 * 基因推荐
 * 
 */

//获取gene sets的信息
export function queryGeneSets(JobId){
    return axios({
        method: "post",
        url: "api/queryGeneSets",
        data: {
            'JobId':JobId,
        },
    })
}

//enricher
export function queryEnricher(JobId,ViewId,organism,gene_set_name,cluster_id,logfoldchanges_threshold,p_threshold,top){
    return axios({
        method: "post",
        url: "api/queryEnricher",
        data: {
            JobId,
            ViewId,
            organism,
            gene_set_name,
            cluster_id,
            logfoldchanges_threshold,
            p_threshold,
            top
        },
    })
}

//gsea
export function queryGsea(JobId,ViewId,organism,gene_set_name,cluster_id,logfoldchanges_threshold,q_threshold,p_threshold,top){
    return axios({
        method: "post",
        url: "api/queryGSEA",
        data: {
            JobId,
            ViewId,
            organism,
            gene_set_name,
            cluster_id,
            logfoldchanges_threshold,
            q_threshold,
            p_threshold,
            top
        },
    })
}


/**
 * socket相关
 */

export function testConnection(){
    return io(window.location.origin,{
        path:window.location.pathname + 'socket.io',
        reconnection: true, // 默认为true
        reconnectionDelay: 100, // 重新连接的延迟（毫秒）
        auth:{
            'type':'Test',
        }

    })
}

export function createNewSocket(JobId){
    return io(window.location.origin,{
        path:window.location.pathname + 'socket.io',
        reconnection: true, // 默认为true
        reconnectionDelay: 100, // 重新连接的延迟（毫秒）
        auth:{
            'type':'Main',
            'JobId':JobId,
        }
    })



}

/**
 * data set相关
 */
export function checkDataSet(JobId,dataset_name){
    return axios({
        method: "post",
        url: "api/check_dataset",
        data: {
            JobId,
            dataset_name,
        },
    })
}



/**
 * 
 * 其他
 * 
 */

// 推荐基因
export function recommendGene(JobId,ViewId,mode){
    return axios({
        method:"post",
        url:"api/recommendGene",
        data:{
            'JobId':JobId,
            'ViewId':ViewId,
            'mode':mode,
        }
    })
}

//删除细胞
export function updateDeleteCells(JobId,ViewId,chosenData){
    return axios({
        method: "post",
        url: "api/updateDeleteCells",
        data: {
            'JobId':JobId,
            'ViewId':ViewId,
            'chosenData':chosenData, 
        },
    });
}

//提交意见
export function sendMessage(form){
    return axios({
        method: "post",
        url: "api/sendMessage",
        data: form,
    });
}


//当前操作实例关闭
export function InstanceClose(){
    
    let reqParams = {};

    return navigator.sendBeacon('/InstanceClose',JSON.stringify(reqParams))
}
