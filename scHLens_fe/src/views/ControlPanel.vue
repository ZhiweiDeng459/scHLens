<template>
    <div>
        <div class="control-panel-container">
            <div class="main-control-panel">
                <basic-control style="height:100%" @startPipeline="startPipeline"></basic-control>
            </div>
        </div>


        <!--加载遮罩层-->
        <LoadingModal ref="PipelineLoading"></LoadingModal>

    </div>
</template>

<script>
import Vue from "vue";
import axios from "axios";
import FileSaver from "file-saver";
import { Form, Tabs, TabPane, Loading, Drawer,Message,Input,Alert,MessageBox} from "element-ui";
import {runPipeline,fetchViewData} from "@/utils/interface"
import BasicControl from "@/components/ControlPanel/BasicControl";
import eventBus from '@/utils/eventBus.js'
import Clipboard from 'v-clipboard'
import LoadingModal from "@/components/LoadingModal";

import io from 'socket.io-client'

Vue.use(Clipboard)

import * as d3 from "d3";
import ElementUI from 'element-ui';
Vue.component(Form.name, Form);
Vue.component(Tabs.name, Tabs);
Vue.component(TabPane.name, TabPane);
Vue.component(Drawer.name, Drawer);
Vue.component(Alert.name, Alert);
Vue.component(MessageBox.name, MessageBox);


Vue.prototype.$message = Message;



export default {
    name: "ControlPanel",
    components: {
        BasicControl,
        LoadingModal,

            },
    props: ["layer"],
    data() {
        return {
            activeName: "Basic",
            dataSet: "",
            mode: "global",
        };
    },
    computed: {
        curData() {
            return this.$store.state.curData;
        },
        dataList() {
            return this.$store.state.dataList;
        },
        JobId(){
            return this.$store.state.JobId
        },
        socketIns(){
            return this.$store.state.socketIns;
        }
    },
    watch: {

    },
    mounted(){
        let self = this
        eventBus.$on('receivedTreeData',(Treedata)=>{//在load job时，改事件会被触发 -> 遍历
            function traverse(root){
                self.createNewViewData(root.data);
                if(root['children'].length != 0){
                    for(let child of root['children']){
                        traverse(child)
                    }
                }
            }
            traverse(Treedata)
            this.$store.commit("toggleCurData", Treedata.ViewId);
        })
    },
    methods: {
        startPipeline(params) {//已经收集到了参数，开始执行pipeline

            console.log('params:',params)

            //监听并实时改进Loading
            this.socketIns.on('get_pipeline_schedule',(data)=>{
                this.$refs['PipelineLoading'].setBottomInfo(`Current Step : ${data['status']}`)
                this.$refs['PipelineLoading'].setPercentage(data['percentage'] * 100)

            })

            //init Loading
            this.$refs['PipelineLoading'].setVisible(true)
            this.$refs['PipelineLoading'].setTopInfo('Loading...')
            this.$refs['PipelineLoading'].setBottomInfo('Preparing...')
            this.$refs['PipelineLoading'].setPercentage(0)


            //排除不合法的局部投影
            if ('local' in params['type'] && params['type']['local']["chosenCells"].length == 0) {
                this.$refs['PipelineLoading'].setVisible(true)
                this.$message({
                    'message':'Please choose cells before createing local plot',
                    'type':'error',
                    'duration':10000,
                    'showClose':true,
                })
                return;
            }

            console.time('single_pipeline')

            runPipeline(params)
                .then((pipeline_response) => {

                    

                    fetchViewData(pipeline_response.data.JobId,pipeline_response.data.ViewId)
                        .then((fetch_response)=>{
                            console.log(fetch_response.data)
                            this.createNewViewData(fetch_response.data)
                            this.$message({
                                'message':'The Pipeline Finished',
                                'type':'success',
                                'showClose':true,
                            })
                            this.$refs['PipelineLoading'].setVisible(false)
                            console.timeEnd('single_pipeline')
                        })
                        .catch((fetch_err)=>{
                            console.log(fetch_err)
                            MessageBox.alert(
                                `<strong>Location：</strong>Unknown
                                    <br>
                                <strong>Details：</strong>A unlocatable error ocurred during fetching view data. Maybe you should check the format of data set or the pipeline configs.`,
                                'Error',
                                {
                                    dangerouslyUseHTMLString: true, 
                                    type: 'error', 
                                }
                            );

                            this.$refs['PipelineLoading'].setVisible(false)
                            console.timeEnd('single_pipeline')
                        })
                })
                .catch((pipeline_err) => {
                    console.timeEnd('single_pipeline')
                    //提示结果
                    let err_data = pipeline_err.response.data;
                    if(err_data.type == 'pipelineException'){
                        MessageBox.alert(
                            `<strong>Location：</strong>${err_data.attach.location}
                                <br>
                            <strong>Details：</strong>${err_data.attach.message}
                                <br>
                            <strong>Advice：</strong>${err_data.attach.advice}`,
                            'Pipeline Error',
                            {
                                dangerouslyUseHTMLString: true,
                                type: 'error',
                            }
                        );

                    }
                    else{
                        MessageBox.alert(
                            `<strong>Location：</strong>Unknown
                                <br>
                            <strong>Details：</strong>A unlocatable error ocurred in running process of the pipeline. Maybe you should check the format of data set or the pipeline configs.`,
                            'Pipeline Error',
                            {
                                dangerouslyUseHTMLString: true, 
                                type: 'error', 
                            }
                        );
                    }
                    this.$refs['PipelineLoading'].setVisible(false)

                });
        },

        createNewViewData(data){//用新增的数据创建视图

            let params = data.paramsObj

            //构建层次化结构
            let node = {}
            if(data.ParentId == 'root'){//根节点
                node = {
                    id: data.ViewId,
                    type: params.type, //目前有global、local
                    layer: 1,
                    Parent:null,
                    children: [],
                };
            }
            else{//非根节点
                //查找父节点
                let ParentId = data.ParentId;
                let ParentData = null;
                for(let d of this.dataList){
                    if(d.ViewId == ParentId){
                        ParentData = d;
                    }
                }
                let ParentNode = ParentData.TreeNode;

                node = {//子节点
                    id: data.ViewId,
                    type: params.type, //目前有global、local
                    layer: ParentNode.layer + 1,
                    Parent:ParentNode,
                    children: [],
                };
            }

            data.TreeNode = node;


            //提交树结构
            this.$store.commit("addNodeToTree",node);

            //提交data数据
            this.$store.commit("addDataObj", data);
            //切换视图
            this.$store.commit("toggleCurData",data.ViewId)
            

            //激活视图（根据实际的流水线）
            //装入activeFlag
            data.activeFlag = {
                'CellProjection':false,
                'GeneProjection':false,
                'GeneExpression':false,
                'MarkerGene':false,
                'TrajectoryInference':false,
                'CellChat':false,
            };
            data.activeFlag['CellProjection'] = true;
            data.activeFlag['GeneProjection'] = true;
            data.activeFlag['GeneExpression'] = true;
            data.activeFlag['MarkerGene'] = true;
            data.activeFlag['TrajectoryInference'] = true;
            data.activeFlag['CellChat'] = true;    

            //清理原有的视图
            this.$store.commit("cleanChooseView")
            //引入默认的散点图
            this.$store.commit("chooseView",[0,"CellProjection"])

            console.log('curData:',this.curData)
        },


    },
};
</script>

<style scoped lang="less">
.control-panel-container{
    position: relative;
    width: 100%;
    height: 100%;
    border-right: 2px solid rgb(200, 200, 200);
    display: flex;
    flex-direction: column;

    .main-control-panel {
        flex: 1 1 0;
        width: 100%;
        /deep/ .el-tabs__item {
            padding: 0 20px;
        }
        .tab {
            margin-top: 10px;
            margin-bottom: 5px;
        }
    }
}

</style>
