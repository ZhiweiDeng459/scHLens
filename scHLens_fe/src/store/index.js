import Vue from "vue";
import Vuex from "vuex";


Vue.use(Vuex);
import { fetchDatasets} from "@/utils/interface";
export default new Vuex.Store({



    state: {
        //protype


        //JobId
        JobId:'',

        dataList: [],
        curData: {
            /**
             * 
             * format of curData:
             *
             * chosenData
             * 
             * activeFlag
             * 
             * ViewId
             * 
             */
            chosenData: [], 
            TreeNode: { layer: 1 } ,
            //视图激活
            activeFlag:{
                'CellProjection':false,
                'GeneProjection':false,
                'GeneExpression':false,
                'MarkerGene':false,
                'TrajectoryInference':false,
                'CellChat':false,
            }},
        curGeneName:[],
        repaintTag:{
            //调用某个视图的部分重绘
            'all':0,
            'Scatter':0,
            'GeneScatter':0,
            'DensityScatter':0,
            'DotPlot':0,
            'Violin':0,
            'ProjectionTree':0,
            'CellChatChord':0,
            'CellChatEdge':0,
        },
        recommendMode:"HighlyVariable",//推荐策略
        pipelineEntityArray:[],//流水线信息

        //数据集相关
        dataSetOptions: [
        ],

        //层次化结构相关
        projectTree: {
            'root':null,
        },
        mergeRoot:undefined,
        chosenViews:[],


        //socket相关
        socketIns:null,//socketIO实例

        //global InfoPanel相关
        infoPanel:null,//infoPanel实例

        //global MessageBoard相关
        messageBoard:null,


        //视图选择相关
        chooseFlags:[{//容器0
                'CellProjection':false,
                'GeneProjection':false,
                'GeneExpression':false,
                'MarkerGene':false,
                'TrajectoryInference':false,
                'CellChat':false,
            },
            {//容器1
                'CellProjection':false,
                'GeneProjection':false,
                'GeneExpression':false,
                'MarkerGene':false,
                'TrajectoryInference':false,
                'CellChat':false,
        }]
    },
    mutations: {
        //设置JobId
        setJobId(state,JobId){
            state.JobId = JobId;
        },
        //添加数据
        addDataObj(state, newData) {
            
            state.dataList.push(newData);            
            //设置默认基因
            state.curGeneName = [newData.defaultGene];
        },
        //删除数据
        deleteDataObj(state,ViewId){
            //如果为删除节点为根节点，那么终止该操作
            if(state.projectTree['root'].id == ViewId){
                return;
            }

            function check_children(root,find_id){//校验是否有继承关系（包括同一节点）
                if (root.id == find_id){
                    return true
                }
                else{
                    for(let child of root['children']){
                        if(check_children(child,find_id))
                            return true
                    }
                }
                return false
            }

            let data = state.dataList.filter((item) => item.ViewId === ViewId)[0];
            let node = data['TreeNode']
            //整理要删除的节点（View及其子集）
            let toDelete = []
            for(let _data of state.dataList){
                if(check_children(node,_data['ViewId']))
                    toDelete.push(_data)
            }
            let toDeleteId = toDelete.map(v=>v.ViewId)

            //从projectionTree中删除
            node['Parent']['children'].splice(node['Parent']['children'].indexOf(node),1)

            //从DataList中删除
            for(let _data of toDelete){
                let index = state.dataList.indexOf(_data)
                state.dataList.splice(index,1)
            }
            //删除先前数据 注意：这里可能出现本该有mergeRoot出现，但是又因为置为了undefined的情况。所以我们在navigation view中有一个resetMergeRoot
            for(let id of toDeleteId){
                if(id == state.mergeRoot){
                    state.mergeRoot = undefined
                    break
                }
            }
            for(let id of toDeleteId){
                let index = state.chosenViews.indexOf(id)
                if(id != -1){
                    state.chosenViews.splice(index,1)
                }
            }
            //重设当前视图
            if(toDelete.includes(state.curData)){
                if(state.dataList.length != 0){
                    //清理前一个curData的相关数据
                    state.curData.chosenData = []
                }
                state.curData = state.dataList.filter((item) => item.ViewId === state.projectTree['root'].id)[0];
                //整理前一个curData的相关数据
                state.curGeneName = [state.curData.defaultGene]        
            }


        },


        //切换数据
        toggleCurData(state, ViewId) {
            if(state.dataList.length != 0){
                //清理前一个curData的相关数据
                state.curData.chosenData = []
            }
            state.curData = state.dataList.filter((item) => item.ViewId === ViewId)[0];
            //整理前一个curData的相关数据
            state.curGeneName = [state.curData.defaultGene]
        },

        //更新选择数据
        updateChosenData(state, chosenData) {
            state.curData.chosenData = chosenData;
        },

        /**curGeneName相关 */
        updateCurGeneName(state, curGeneName) {//直接修改
            state.curGeneName = curGeneName;
        },
        addToCurGeneName(state, geneName) {//添加基因
            if(!state.curGeneName.includes(geneName)){
                state.curGeneName.unshift(geneName)
            }
        },
        addGroupToCurGeneName(state, geneNames) {//添加基因一组
            geneNames.forEach(geneName=>{
                if(!state.curGeneName.includes(geneName)){
                    state.curGeneName.unshift(geneName)
                }
            })
        },
        deleteFromCurGeneName(state,geneName){//删除指定基因
            if(state.curGeneName.includes(geneName) && state.curGeneName.length >1){
                state.curGeneName.splice(state.curGeneName.indexOf(geneName),1)
            }
        },
        initCurGeneName(state){//恢复curGeneName
            state.curGeneName = [state.curData.defaultGene]
        },


        //通过新增节点调整树
        addNodeToTree(state,node) {
            if(node.layer == 1){//根节点-新建树
                state.projectTree['root'] = node;

            }
            else{//子节点-把子节点添加到树的对应父节点的children下
                //找寻父节点
                state.dataList.find(d=>{
                    return d.ViewId == node.Parent.id
                }).TreeNode.children.push(node);

            }
        },
        
        //更新mergeRoot
        updateMergeRoot(state,mergeRoot){
            state.mergeRoot = mergeRoot;
        },

        //更新视图的新数据
        updateViewData(state,data){ 
            let ViewId = data.ViewId;
            let dataListIndex = state.dataList.map(d=>d.ViewId).indexOf(ViewId);

            //继承数据
            data.TreeNode = state.dataList[dataListIndex].TreeNode
            data.activeFlag = state.dataList[dataListIndex].activeFlag

            if(state.curData === state.dataList[dataListIndex]){//如果当前数据即是要替换的数据
                state.dataList[dataListIndex] = data;
                state.curData = data
                state.curData.chosenData = []
                state.curGeneName = [state.curData.defaultGene]
            }
            else{
                state.dataList[dataListIndex] = data;
            }
        },

        /**
         * 
         * 视图管理（激活、关闭）
         *
         */
        //为图框选择视图
        chooseView(state,params){
            state.chooseFlags[params[0]][params[1]] = true
        },
        //取消图框中对某个视图的选择
        unchooseView(state,params){
            state.chooseFlags[params[0]][params[1]] = false;
        },
        //取消所有图框中视图的选择
        cleanChooseView(state){
            state.chooseFlags.forEach((chooseFlag)=>{
                for(let key in chooseFlag){
                    chooseFlag[key] = false;
                }
            })
        },
        
        /**
         * 
         * 更新数据集
         * 
         */
        updateDatasets(state,JobId){
            fetchDatasets(JobId).then((response)=>{
                let datasetConfigs = response.data
                
                //整理格式
                let dataSetOptions = []
                for(let i = 0;i < datasetConfigs.length;i++){
                    let option = {
                        index:i,
                        value:datasetConfigs[i],
                        label:datasetConfigs[i].name + (datasetConfigs[i].from=='sample'?' (sample)':'')
                    }
                    dataSetOptions.push(option)
                }
                state.dataSetOptions = dataSetOptions;
            }).catch((err) => {
                console.log(err)
            })
        },



        /**
         * socket相关
         */
        setSocket(state,socketIns){
            state.socketIns = socketIns
        },
        clearSocket(state,socketIns){
            if(state.socketIns !== null)
                state.socketIns.disconnect();
            state.socketIns = null;
        },

        /**
         * infoPanel相关（游走小信息板）
         */
        setInfoPanel(state,infoPanel){
            state.infoPanel = infoPanel
        },

        /**
         * MessageBoard相关（固定大信息板）
         */
        setMessageBoard(state,messageBoard){
            state.messageBoard = messageBoard
        },
        /**
         * AnnoRecomAddiInfoPanel相关（标注额外信息板）
         */
        setAnnoRecomAddiInfoPanel(state,AnnoRecomAddiInfoPanel){
            state.AnnoRecomAddiInfoPanel = AnnoRecomAddiInfoPanel
        },

        /**
         * 其他
         */
        //刷新视图
        refreshView(state,type){
            if(type == 'all'){
                for(let key of Object.keys(state.repaintTag)){
                    state.repaintTag[key] = state.repaintTag[key] + 1;
                }
            }
            else{
                state.repaintTag[type] = state.repaintTag[type] + 1;
            }
            
        }

    },
    actions: {},
    modules: {},
});
