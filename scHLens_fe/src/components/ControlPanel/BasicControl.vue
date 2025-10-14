<template>
    <div>
        <el-form class="form">
            <b style="color:white;font-size:20px">Job Id</b>
            <HR style="margin:1px 0px"></HR>
            <el-form-item class="form-item" style="margin:5px 0px 10px 0px;">
                <div
                    style="
                    display:flex;
                    align-items:center;
                    justify-content:space-between;
                    width:100%">
                    <el-input
                        v-model="JobId"
                        size="small"
                        style="flex:1 1 0;margin-right:10px"
                        :disabled="true"
                    >
                    <el-button
                        type="primary"
                        icon="el-icon-document-copy"
                        slot="append"
                        @click="copyJobId"
                        ></el-button>
                    </el-input>

                </div>
            </el-form-item>

            <b style="color:white;font-size:20px">Data Set</b>
            <HR style="margin:1px 0px"></HR>
            <el-form-item class="dataset-item" style="margin-bottom:0px">
                <div
                    style="
                    display:flex;
                    align-items:center;
                    justify-content:space-between;
                    width:100%">
                    <el-select
                        v-model="dataset"
                        placeholder="Select..."
                        size="small"
                        style="flex:1 1 0;margin-right:10px"
                        >
                        <el-option v-for="item in dataSetOptions" :key="item.index" :label="item.label" :value="item">
                            <div style="display: flex;flex-direction: row;justify-content: space-between;align-items: center;">
                                <span>{{item.label}}</span>
                            </div>
                        </el-option>
                    </el-select>
                    <el-button
                        type="primary"
                        icon="el-icon-plus"
                        size="medium"
                        style="width:26px;height:26px;padding:0px;margin-right:10px;"
                        @click="openUploadPanel"></el-button>
                </div>
                <div style="width:100%;line-height: 32px;">
                    <el-input v-model="curDataSetCellNum" :readonly="true" size="mini" style="width: 95%;">
                        <div slot="append" style="width:40px;display: flex;align-items: center;justify-content: center;">
                            <a style="font-size:12px;font-family:YaHei;">Cells</a>
                        </div>
                    </el-input>
                </div>
                <div style="width:100%;line-height: 32px;">
                    <el-input v-model="curDataSetGeneNum" :readonly="true" size="mini" style="width: 95%;">
                        <div slot="append" style="width:40px;display: flex;align-items: center;justify-content: center;">
                            <a style="font-size:12px;">Genes</a>
                        </div>
                    </el-input>
                </div>
            </el-form-item>


            <el-form-item class="button-list">
                <el-button type="primary" class="form-buttons" size="mini" @click="callGlobalPipelineConfig">Set Pipeline</el-button>
                <el-popover
                    placement="bottom"
                    width="300"
                    trigger="click">
                    <GeneFilter/>
                    <el-button slot="reference" type="primary" class="form-buttons" size="mini">Filter Cells</el-button>
                </el-popover>
                <el-button type="primary" class="form-buttons" size="mini" @click="exportJob">Export Job</el-button>
                <el-popconfirm
                    confirm-button-text='Yes'
                    cancel-button-text='No'
                    icon="el-icon-info"
                    icon-color="red"
                    confirm-button-type="Primary"
                    cancel-button-type="Primary"
                    @confirm="deleteJob"
                    title="Are you sure you want to delete the current job? This action is irreversible.">
                    <el-button type="danger" class="form-buttons" size="mini" slot="reference">Delete Job</el-button>
                </el-popconfirm>

            </el-form-item>

            <div style="flex:1 1 0"></div>
            
            <el-form-item class="end-list">
                <el-button type="danger" class="form-buttons" size="mini" @click="ExitJob">Exit Job</el-button>
            </el-form-item>
        </el-form>

        <!--pipeline对话框-->
        <el-dialog
            title="The Pipeline"
            :visible.sync="pipelineDialogVisible"
            :close-on-click-modal="false"
            :modal-append-to-body="false"
            width="90%"
            top="5vh">
            <div style="display:flex"> 
                <el-popover
                    placement="right"
                    width="400px"
                    trigger="click"
                    v-model = "pipelineAddPopVisble">
                    <div v-for="Entity in pipelineEntity" 
                        :key="Entity.id" 
                        :class="{'ToAddEntity':!InPipelineEntityArray(Entity.id),'AddedEntity':InPipelineEntityArray(Entity.id)}" 
                        @mouseenter="enterAddEntity" 
                        @mouseleave="leaveAddEntity" 
                        @click="AddEntity(Entity)" 
                        style="width:300px;margin:10px 5px; display:flex; justify-content:space-between; align-items:center"
                        >
                        <b style="color:black">{{Entity.name}}</b>
                        <i class="el-icon-check" v-show="InPipelineEntityArray(Entity.id)"></i>
                    </div>
                    <el-button slot="reference" icon="el-icon-plus" style="width:70px;height:30px;padding:0px">Add</el-button>
                </el-popover>
                <el-button icon="el-icon-plus" style="width:70px;height:30px;padding:0px;margin:0px 5px" @click="autoConfig">All</el-button>
                <el-button icon="el-icon-delete" style="width:70px;height:30px;padding:0px;margin:0px 5px" @click="clearEntity">Clear</el-button>
                <div style="flex:1 1 0"></div>
                <el-button icon="el-icon-check" type="success" style="width:100px;height:30px;padding:0px;margin:0px 5px" @click="startPipeline()">Run</el-button>
            </div>
            <div style="display:flex;flex-direction:row;overflow-X:auto">
                <el-container v-for="element in pipelineEntityArray" :key="element.index" class="StepContainerInPipeline" style="flex-direction:column">
                    <el-header height="30px" style="padding:0px 4px;border-radius:4px;display:flex;flex-direction:row;justify-content:space-between">
                            <b style="font-size:20px;margin-top:5px;display:block">{{element.name}}</b>
                            <el-button type="danger" circle icon="el-icon-close" style="padding:3px;margin:3px 0px" @click="deleteEntity(element)"></el-button>
                    </el-header>
                    <el-main style="padding:0px 0px">
                        <SamplingParams v-if="element.id == 'SP'" :ref="'PipelineEntity' + element.index"/>
                        <QualityControlParams v-if="element.id == 'QC'" :ref="'PipelineEntity' + element.index"/>
                        <TransformationParams v-else-if="element.id == 'TS'" :ref="'PipelineEntity' + element.index" :mode="mode"/>
                        <GeneSelectionParams v-else-if="element.id == 'FS'" :ref="'PipelineEntity' + element.index"/>
                        <DimensionReductionParams v-else-if="element.id == 'DR'" :ref="'PipelineEntity' + element.index"/>
                        <ClusteringParams v-else-if="element.id == 'CL'" :ref="'PipelineEntity' + element.index"/>
                        <MarkerGenesParams v-else-if="element.id == 'MK'" :ref="'PipelineEntity' + element.index"/>
                        <!-- <TrajectoryInferenceParams v-else-if="element.id == 'TI'" :ref="'PipelineEntity' + element.index"/>
                        <DataIntegrationParams v-else-if="element.id == 'DI'" :ref="'PipelineEntity' + element.index" :ref_dataset="dataset"/> -->
                        <CellChatParams v-else-if="element.id == 'CC'" :ref="'PipelineEntity' + element.index"/>
                    </el-main>
                </el-container>
            </div>
        </el-dialog>

        <!--上传数据集对话框-->
        <el-dialog
            title="Upload the Dataset"
            :visible.sync="uploadDialogVisible"
            :close-on-click-modal="false"
            width="400px">
            <div style="
                display:flex;
                width:100%;
                height:100%;
                flex-direction:column;
                align-items:center
            ">
                <div style="
                     display:flex;
                     flex-direction:column;
                     width: 360px;
                     margin: 10px 0;
                     justify-content:space-between;
                     background-color:#fff;
                     border-radius:6px;
                     border:1px solid #eaeaea;
                     ">
                    <div style="text-align:center;margin-top:10px;">
                        <a style="font-size:17px;font-family:YaHei">Requirements</a>
                    </div>

                    <div class="upload-require-container" v-if="upType=='csv'">
                        <el-steps direction="vertical">
                            <el-step v-for="item in fileRequirements['csv']" :key="item['label']" :title="item['label']" :status="item['exist']?'success':'wait'"></el-step>
                        </el-steps>
                    </div>
                    <div class="upload-require-container" v-if="upType=='10x-mtx'">
                        <el-steps direction="vertical">
                            <el-step v-for="item in fileRequirements['10x-mtx']" :key="item['label']" :title="item['label']" :status="item['exist']?'success':'wait'"></el-step>
                        </el-steps>
                    </div>
                    <div class="upload-require-container" v-if="upType=='h5ad'">
                        <el-steps direction="vertical">
                            <el-step v-for="item in fileRequirements['h5ad']" :key="item['label']" :title="item['label']" :status="item['exist']?'success':'wait'"></el-step>
                        </el-steps>
                    </div>
                    <div class="upload-require-container" v-if="upType=='h5df'">
                        <el-steps direction="vertical">
                            <el-step v-for="item in fileRequirements['h5df']" :key="item['label']" :title="item['label']" :status="item['exist']?'success':'wait'"></el-step>
                        </el-steps>
                    </div>
                </div>

                <el-upload
                    drag
                    multiple
                    ref="upload"
                    :file-list="this.fileList"
                    :limit="fileRequirements[upType].length"
                    :on-change="changeUploadFile"
                    :on-remove="removeUploadFile"
                    :on-success="successUploadFile"
                    :on-error="errorUploadFile"
                    :on-exceed="exceedUploadFile"
                    :data="{
                        'name':uploadDatasetName,
                        'type':upType,
                        'JobId':JobId
                    }"
                    :auto-upload="false"
                    :action="uploadPath">
                    <i class="el-icon-upload"></i>
                    <div class="el-upload__text">Drag file to this,or &nbsp;<em>Click to add</em></div>
                </el-upload>


                <div class="el-upload__tip" slot="tip" style="font-family:YaHei">Set dataset type and name:</div>
                <div
                    style="
                    display:flex;
                    width: 360px;
                    margin-top:10px;
                    justify-content:space-between
                    ">

                    <el-input
                        placeholder="Name"
                        v-model="uploadDatasetName"
                        style="width:220px">
                        <el-select
                                v-model="upType"
                                placeholder="Select..."
                                size="small"
                                style="width:100px;"
                                slot="prepend"
                                >
                                <el-option v-for="item in dataTypeOptions" :key="item.index" :label="item.label" :value="item.value"> </el-option>
                        </el-select>
                    </el-input>
                        <el-button size="small" type="success" @click="submitUpload">upload to server</el-button>
                </div>
            </div>
        </el-dialog>

        
    </div>
</template>

<script>
import Vue from "vue";
import {DescriptionsItem,Descriptions,Popconfirm, Upload,Form, FormItem, Button, Select, Option, Radio, Tooltip, Dialog, Container, Main, Aside, Header,Footer,Popover,Loading,MessageBox} from "element-ui";

import QualityControlParams from "@/components/ControlPanel/PipelineEntity/QualityControlParams";
import TransformationParams from "@/components/ControlPanel/PipelineEntity/TransformationParams";
import SamplingParams from "@/components/ControlPanel/PipelineEntity/SamplingParams"
import DimensionReductionParams from "@/components/ControlPanel/PipelineEntity/DimensionReductionParams";
import ClusteringParams from "@/components/ControlPanel/PipelineEntity/ClusteringParams";
import MarkerGenesParams from "@/components/ControlPanel/PipelineEntity/MarkerGenesParams";
import GeneSelectionParams from "@/components/ControlPanel/PipelineEntity/GeneSelectionParams";
import TrajectoryInferenceParams from "@/components/ControlPanel/PipelineEntity/TrajectoryInferenceParams"
import DataIntegrationParams from "@/components/ControlPanel/PipelineEntity/DataIntegrationParams"
import CellChatParams from "@/components/ControlPanel/PipelineEntity/CellChatParams"
import GeneFilter from "@/components/GeneFilter"
import eventBus from '@/utils/eventBus.js'
import axios from 'axios';
import { saveAs } from 'file-saver';

import {exportJob,checkDataSet,clearCallback,deleteJob} from "@/utils/interface.js"


import Clipboard from 'v-clipboard'
Vue.use(Clipboard)

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Button.name, Button);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Dialog.name, Dialog);
Vue.component(Popover.name, Popover);

Vue.component(Container.name,Container);
Vue.component(Main.name,Main);
Vue.component(Aside.name,Aside);
Vue.component(Header.name, Header);
Vue.component(Footer.name, Footer);
Vue.component(Upload.name,Upload);
Vue.component(MessageBox.name,MessageBox);
Vue.component(Popconfirm.name,Popconfirm);

Vue.component(Descriptions.name,Descriptions)
Vue.component(DescriptionsItem.name,DescriptionsItem)
export default {
    name: "BasicControl",

    components:{
        SamplingParams,
        QualityControlParams,
        ClusteringParams,
        DimensionReductionParams,
        TransformationParams,
        MarkerGenesParams,
        GeneSelectionParams,
        // TrajectoryInferenceParams,
        // DataIntegrationParams,
        CellChatParams,
        GeneFilter,
    },

    data() {
        return {
            /**
             * 数据集相关
             */
            dataset: {},
            mode: "global",
            isSelect: false,
            curDataSetCellNum:0,
            curDataSetGeneNum:0,


            /**
             * 文件上传相关
             */
            uploadPath : '/api/upload',
            uploadDialogVisible : false,
            fileList:[],
            upType:'csv',
            uploadDatasetName:'',
            dataTypeOptions:[
                {
                    index: 1,
                    value: "csv",
                    label: "csv",
                },
                {
                    index: 2,
                    value: "10x-mtx",
                    label: "10x-mtx",
                },
                {
                    index: 3,
                    value: "h5ad",
                    label: "h5ad",
                },
                {
                    index: 4,
                    value: "h5df",
                    label: "h5df",
                },
            ],
            fileRequirements:{/**每种类型的数据集需求的文件信息 */
                'csv':[
                    {
                        'label':'*.csv',
                        'check_format':/^.*\.csv$/i,//check_format是指用来匹配上传文件名是否满足的正则字面量
                        'exist':false,
                    }
                ],
                '10x-mtx':[
                    {
                        'label':'genes.tsv',
                        'check_format':/^genes\.tsv$/i,
                        'exist':false,

                    },{
                        'label':'barcodes.tsv',
                        'check_format':/^barcodes\.tsv$/i,
                        'exist':false,
                    },{
                        'label':'matrix.mtx',
                        'check_format':/^matrix\.mtx$/i,
                        'exist':false,
                    }
                ],
                'h5ad':[
                    {
                        'label':'*.h5ad',
                        'check_format':/^.*\.h5ad$/i,
                        'exist':false,
                    }
                ],
                'h5df':[
                    {
                        'label':'*.h5',
                        'check_format':/^.*\.h5$/i,
                        'exist':false,
                    }
                ],
                
            },


            /**
             * pipline相关
             */
            pipelineDialogVisible : false,
            pipelineAddPopVisble : false,
            pipelineRawEntity:[ //保存所有的的pipeline entity，
                {
                    rank:0,
                    id:'SP',
                    name:"Downsampling"
                },
                {
                    rank:1,
                    id:'QC',
                    name:"Quality Control"
                },
                {
                    rank:2,
                    id:'TS',
                    name:"Normalization", 
                },
                {
                    rank:3,
                    id:'FS',
                    name:"Gene Selection", 
                },
                // {
                //     rank:2,
                //     id:'DI',
                //     name:"Data Integration"
                // },
                // {
                //     rank:3,
                //     id:'NB',
                //     name:'Neighbor'
                // },
                {
                    rank:4,
                    id:'DR',
                    name:'Visualization',
                },
                {
                    rank:5,
                    id:'CL',
                    name:'Clustering',
                },
                {
                    rank:6,
                    id:'MK',
                    name:"DEG Identification"
                },
                // {
                //     rank:7,
                //     id:'TI',
                //     name:"Trajectory Inference"
                // },
                {
                    rank:7,
                    id:'CC',
                    name:"Cell Communication"
                }
                ],
            pipelineEntity: [ //保存当前状态下能使用的pipelineEntity
            {
                    rank:0,
                    id:'SP',
                    name:"Downsampling"
                },
                {
                    rank:1,
                    id:'QC',
                    name:"Quality Control"
                },
                {
                    rank:2,
                    id:'TS',
                    name:"Normalization", 
                },
                {
                    rank:3,
                    id:'FS',
                    name:"Gene Selection", 
                },
                // {
                //     rank:2,
                //     id:'DI',
                //     name:"Data Integration"
                // },
                // {
                //     rank:3,
                //     id:'NB',
                //     name:'Neighbor'
                // },
                {
                    rank:4,
                    id:'DR',
                    name:'Visualization',
                },
                {
                    rank:5,
                    id:'CL',
                    name:'Clustering',
                },
                {
                    rank:6,
                    id:'MK',
                    name:"DEG Identification"
                },
                // {
                //     rank:7,
                //     id:'TI',
                //     name:"Trajectory Inference"
                // },
                {
                    rank:7,
                    id:'CC',
                    name:"Cell Communication"
                }
                ],
    
            /**
             * 基因推荐相关
             */
            recommendMode: "HighlyVariable",
            recommendModeOptions: [
                {
                    index: 0,
                    label: "HighlyVariable",
                },
                {
                    index: 1,
                    label: "CellMarker",
                },
            ],


        };
    },
    computed: {
        curData() {
            return this.$store.state.curData;
        },
        chosenData() {
            return this.$store.state.curData.chosenData;
        },
        pipelineEntityArray(){ //当前用户已经使用的params entity
            return this.$store.state.pipelineEntityArray;
        },
        dataSetOptions(){
            return this.$store.state.dataSetOptions;
        },
        JobId(){
            return this.$store.state.JobId
        },

    },
    watch: {
        chosenData() {
            //当有数据被选中，那么数据集不能再被切换
            if (this.chosenData === undefined || this.chosenData === null) {
                this.isSelect = false;
                return;
            }
            if (this.chosenData.length == 0) {
                this.isSelect = false;
            } else {
                this.isSelect = true;
                let datasetParam = this.curData.paramsObj.dataset;
                for (let item of this.dataSetOptions) {
                    if (JSON.stringify(item.value) == JSON.stringify(datasetParam)) {
                        this.dataset = item;
                        break;
                    }
                }
            }
        },
        mode(newValue){
            //当mode切换的时候，pipline的设置受限，仅限与FS以及后续的下游分析
            this.clearEntity()
            if(newValue == 'local'){
                // let access = ['TS','NB','DR','CL','TI','MK','CC']
                let access = ['QC','TS','FS','DR','CL','MK','CC'] //access表示local能包含哪些步骤
                this.pipelineEntity = this.pipelineRawEntity.filter((v)=>{
                    return access.indexOf(v['id']) != -1
                })
            }
            else{
                this.pipelineEntity = [...this.pipelineRawEntity]
            }
        },
        dataset(newValue){//当前选择的数据集转换时
            //切换数据集的size信息
            this.curDataSetCellNum = Number(newValue.value['cell_num'])
            this.curDataSetGeneNum = Number(newValue.value['gene_num'])
            if(isNaN(this.curDataSetCellNum))
                this.curDataSetCellNum = 0
            if(isNaN(this.curDataSetGeneNum))
                this.curDataSetGeneNum = 0
        }
    },
    methods: {
        callGlobalPipelineConfig(){
            eventBus.$emit('callPipelineConfig','global');
        },
        startPipeline() { //主要是前期的参数整理
            //异常处理
            if(Object.keys(this.dataset).length == 0 && this.mode == 'global'){//没有数据集的一场处理
                this.pipelineDialogVisible = false;
                this.$message({
                    'message':'Please choose a data set',
                    'type':'error',
                    'duration':8000,
                    'showClose':true,
                }) 
                return ;
            }
            

            //整理参数
            let Params = {}
            //装入JobId
            Params['JobId'] = this.JobId;
            //装入数据集参数
            Params['dataset'] = this.dataset.value;
            //装入投影种类
            Params['type'] = {}
            Params['type'][this.mode] = {}
            if(this.mode == 'local'){
                Params['type']['local']['chosenCells'] = this.curData.chosenData;
                Params['ParentId'] = this.curData.ViewId;
            }
            else{
                Params['ParentId'] = 'root';
            }
            // //提取各个模块的参数
            // this.pipelineEntityArray.forEach(entity => {
            //     let paramsResult = this.$refs['PipelineEntity' + entity.index][0].getParams();

            //     //参数合法性检查与报错
            //     if('error' in paramsResult){
            //         MessageBox.alert(
            //                 `<strong>Location：</strong>${paramsResult.location}
            //                     <br>
            //                 <strong>Details：</strong>${paramsResult.message}`,
            //                 'Parameter Error',
            //                 {
            //                     dangerouslyUseHTMLString: true,
            //                     type: 'error',
            //                 }
            //             );
            //         return;
            //     }


            //     Params[entity.id] =  paramsResult;    
            // });

            
            for(let entity of this.pipelineEntityArray){
                let paramsResult = this.$refs['PipelineEntity' + entity.index][0].getParams();

                //参数合法性检查与报错
                if('error' in paramsResult){
                    MessageBox.alert(
                            `<strong>Location：</strong>${paramsResult.location}
                                <br>
                            <strong>Details：</strong>${paramsResult.message}`,
                            'Parameter Error',
                            {
                                dangerouslyUseHTMLString: true,
                                type: 'error',
                            }
                        );
                    return;
                }


                Params[entity.id] =  paramsResult;    

            }


            //开始计算
            this.$emit("startPipeline",Params); //提交给上一层进行执行相关的工作

            //关闭pipeline对话框
            this.pipelineDialogVisible = false;
            this.pipelineAddPopVisble = false;

        },
        exportJob() {
            const loading = Loading.service({ fullscreen: true,text:this.loadingText });
            loading.text = `Compressing, please wait...`
            exportJob(this.JobId,(processEvent)=>{
                console.log('processEvent',processEvent)
                // this.loadingText = `已加载： ${processEvent.loaded / processEvent.total}`
                loading.text = `Downloading: ${(processEvent.loaded / processEvent.total * 100).toFixed(3)}%`
            })
                .then((response)=>{
                    const blob = new Blob([response.data])            
                    // const link = document.createElement('a')
                    // link.download = this.JobId + '.zip' // a标签添加属性
                    // link.style.display = 'none'
                    // link.href = URL.createObjectURL(blob)
                    // document.body.appendChild(link)
                    // link.click() // 执行下载
                    // URL.revokeObjectURL(link.href)  // 释放 bolb 对象
                    // document.body.removeChild(link) // 下载完成移除元素
                    saveAs(blob, this.JobId + '.zip');
                    loading.close()
                    clearCallback(this.JobId,'api/exportJob').then(res=>console.log(res)).catch(err=>console.log(err))

                })
                .catch((err)=>{
                    console.log(err)
                    loading.close()
                    clearCallback(this.JobId,'api/exportJob').then(res=>console.log(res)).catch(err=>console.log(err))

                })

        },
        deleteJob(){//删除当前job
            deleteJob(this.JobId).then((res)=>{
                alert("Deleted successfully. The page will refresh.");
                window.location.reload();
            }).catch((err)=>{
                alert("Deletion failed. Please try again.");
            })


        },
        ExitJob(){//退出当前Job
            window.location.reload();
        },
        /**
         * 
         * 流水线相关
         * 
         */

        openPipeline(){//打开自定义流水线窗口
            this.pipelineDialogVisible = true;        
        },
        generatePipelineEntityIndex(){
            //生成流水线实体的Index（随机）
            let length = 6;
            let characterList = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'];
            let existViewIDs = this.pipelineEntityArray.map(item => item['id'])
            let ID = '';
            do{
                ID = ''
                for(let i = 0;i < length; i++){
                    ID += characterList[Math.floor(Math.random()*characterList.length)];
                }
            }
            while(ID in existViewIDs)
            return ID;
        },
        InPipelineEntityArray(id){
            //检查某个id的entity是否在pipelineArray内
            return this.pipelineEntityArray.find(item=>id==item.id)!==undefined
        },
        AddEntity(entityData){
            //当鼠标点击了增添选项 值得注意的是，这里需要排序
            if(this.InPipelineEntityArray(entityData.id))
                return ;
            let self = this;
            this.pipelineAddPopVisble = false;
            this.pipelineEntityArray.push({
                'index':self.generatePipelineEntityIndex(),
                'name':entityData.name,
                'id':entityData.id,
                'rank':entityData.rank,
            });
            this.pipelineEntityArray.sort((a,b)=>a.rank-b.rank)
        },
        AddAllEntity(){
            this.pipelineEntity.forEach(item=>{
                this.AddEntity(item)
            })
        },
        enterAddEntity(event){
            //当鼠标进入在增添选项上时
            if(event.currentTarget.className == 'AddedEntity')
                return;
            event.currentTarget.style.backgroundColor = '#F5F5F5'
            event.currentTarget.childNodes[0].style.color = '#3366cc';
            
        },
        leaveAddEntity(event){
            //当鼠标离开在增添选项上时
            event.currentTarget.style.backgroundColor = 'white'
            event.currentTarget.childNodes[0].style.color = 'black';
        },
        deleteEntity(entityData){
            //删除流水线实体
            this.pipelineEntityArray.splice(this.pipelineEntityArray.indexOf(entityData),1);
        },
        clearEntity(){
            //清空流水线
            this.pipelineEntityArray.splice(0,this.pipelineEntityArray.length)
        },
        autoConfig(){
            //自动设置
            this.clearEntity();
            this.AddAllEntity();
        },
 
        /**
         * 
         * 文件上传相关
         * 
         */
        openUploadPanel(){
           this.uploadDialogVisible = true;
        },
        closeUploadPanel(){
            this.uploadDialogVisible = false;
        },
        submitUpload(){/**文件上传函数 */
            //check1：检查文件需求是否被满足
            for(let f of this.fileRequirements[this.upType]){
                if(!f['exist']){
                    this.$message({
                    'message':'Please select valid files to upload',
                    'type':'error',
                    'showClose':true,
                    })
                    return;              
                }
            }       
            //check2：检查数据集名词
            if(this.uploadDatasetName == ''){
                this.$message({
                    'message':'Please input a valid dataset\'s name',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            //check3：重名检查
            for(let dataset of this.dataSetOptions){
                if(this.uploadDatasetName == dataset.label){
                    this.$message({
                        'message':'This dataset name has been used',
                        'type':'error',
                        'showClose':true,
                    })
                return;
                }
            }
            //check4：文件名非法字符检查
            let forbiddenCharas = ['/','+','-','.']
            for(let c of forbiddenCharas){
                if(this.uploadDatasetName.includes(c)){
                    this.$message({
                        'message':`Invalid name: contains forbidden character "${c}"`,
                        'type':'error',
                        'showClose':true,
                    })
                    return;
                }
            }

            this.$refs.upload.submit();
            this.updateFileRequirements()

        },
        changeUploadFile(file,fileList){
            this.fileList = fileList//由于fileList不会自动更新，所以需要手动更新
            this.updateFileRequirements()
        },
        removeUploadFile(file,fileList){
            this.fileList = fileList//由于fileList不会自动更新，所以需要手动更新
            this.updateFileRequirements()

        },
        successUploadFile(response,file,fileList){ //注意，这里的fileList是上传前的fileList，而不是上传后清空的fileList

            const loading = Loading.service({ fullscreen: true });


            //检查是否所有的文件都被上传
            let allSuccess = true;
            for(let f of fileList){
                if(f.status != 'success'){
                    allSuccess = false;
                    break;
                }
            }
            if(allSuccess){
                checkDataSet(this.JobId,this.uploadDatasetName)
                .then((response)=>{
                    this.$message({
                        'message':'UPLOAD COMPLETE',
                        'type':'success',
                        'showClose':true,
                    })
                    this.$store.commit("updateDatasets",this.JobId)
                    loading.close()

                }).catch((err)=>{
                    this.$message({
                        'message':'CONNECTION LOST OR THIS DATASET IS INVALID',
                        'type':'error',
                        'showClose':true,
                    })
                    this.$store.commit("updateDatasets",this.JobId)
                    loading.close()

                })

                this.$refs.upload.clearFiles();
                this.upType = 'csv'
                this.uploadDatasetName = ''
                this.closeUploadPanel();

                this.fileList = []
                this.updateFileRequirements()

            }

        },
        errorUploadFile(err,file,fileList){ //注意，这里的fileList是上传前的fileList，而不是上传后清空的fileList
            this.$message({
                'message':'UPLOAD FAILED',
                'type':'error',
                'showClose':true,
            })
            this.$refs.upload.clearFiles();
            this.fileList = []
            this.upType = 'csv'
            this.uploadDatasetName = ''
            this.updateFileRequirements()
        },
        exceedUploadFile(files, fileList){//文件超出个数限制
            this.$message({
                'message':`The number of files uploaded cannot exceed the requirement.`,
                'type':'info',
                'showClose':true,
            })
        },
        updateFileRequirements(){//更新上传文件的满足情况

            for(let r_key in this.fileRequirements){
                for(let i = 0;i < this.fileRequirements[r_key].length;i++){
                    this.fileRequirements[r_key][i]['exist'] = false; //防止当fileList为空时，无法将原本的true exist置为false
                    for(let f of this.fileList){
                        if(this.fileRequirements[r_key][i]['check_format'].test(f['name'])){
                            this.fileRequirements[r_key][i]['exist'] = true;
                            break;
                        }
                    }
                }
            }
        },



        /**
         * 
         * 剪贴板相关
         * 
         */
        copyJobId(){//复制Job Id到剪切板
            this.$clipboard(this.JobId);
            this.$message({
                'message':'Copy the Job Id to clipboard successfully',
                'type':'success',
                'showClose':true,
            })
        },

    
    },
    mounted(){
        //绑定eventbus
        eventBus.$on('callPipelineConfig',(mode)=>{//开始计算
            this.mode = mode;//global or local
            this.openPipeline();
        })
        //设定上传数据集文件的url
        this.uploadPath = window.location.pathname + 'api/upload'
    }
};
</script>

<style scoped lang="less">
.form {
    padding: 0px 0 0 15px;
    display: flex;
    flex-direction: column;
    height: 100%;
    .form-item {
        display: flex;
        margin-bottom: 0;
    }
    .dataset-item{
        display: flex;
        flex-direction: column;
    }
    .button-list {
        display: flex;
        flex-direction: column;
        align-items: center;
        .form-buttons {
            width: 91%;
            height: 35px;
            margin: 15px 0 0 0px;
            /deep/ span{
                font-size:16px;
            }
        }
    }
    .end-list {
        display: flex;
        flex-direction: column;
        align-items: center;
        /deep/ .el-form-item__content{
            width: 100%;
        }
        .form-buttons {
            width: 91%;
            height: 35px;
            margin: 15px 0 0 0px;
            /deep/ span{
                font-size:16px;
            }
        }
    }

}

.upload-require-container{
    padding:15px;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}

.ToAddEntity{//添加界面出现的流水线实体
    padding:10px 0px;
    cursor:pointer;
}
.AddedEntity{
    padding:10px 0px;

}

/deep/ .el-dialog__body{//流水线对话框的主体（除去标题）
    background-color:#F5F5F5;
    border-top: 2px solid lightgray;
    padding: 10px 10px;
    

    .StepContainerInPipeline{
    height: 680px;
    flex:0 0 300px;
    margin:10px;
    padding:7px;

    border:2px solid lightgray;
    border-radius: 10px;
    
    background-color :white;
}
}

/deep/  .el-radio-button__inner{ //radio按钮的外界尺寸
    width:110px;
    font-family:YaHei;
}

/deep/ .el-input__inner{

}

/deep/ .el-radio-button .el-radio-button__inner{//radio选择框中字体的样式
    font-size:16px;
    font-family:YaHei;
    padding: 7px 7px;
}

/deep/ .el-select-dropdown__item span{ //选择下拉框中字体的样式
    font-size:16px
}

/deep/ .el-radio-button__inner{//radio选择框中
    width:117px;
    height: 35px;
}

/deep/ .el-upload__text{ //上传文件的中心文本
    font-size:18px;

}
/deep/ .el-upload__tip{ //上传文件的提示文本
    font-size:16px;
}
/deep/ .el-dialog__header span{
    font-size:20px;
}

/deep/ .el-button--small{
    padding:2px 5px;
    span{
        font-size: 16px;
    }
}

/deep/ .el-icon-document-copy{
    color:black;

}


/deep/ .el-input__inner{
    color:#606266 !important;
}
</style>
