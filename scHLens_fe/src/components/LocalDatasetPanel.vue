<template>
        <el-dialog
            :visible.sync="showDialog"
            width="525px"
            :show-close="true"
            :lock-scroll="true"
            :modal="true"
            :close-on-click-modal="false"
            title="Save the Chosen Cells as Dataset"
            >
            <div class="localDatasetContainer">
                <a style="
                    font-size:22px;
                    margin-top:50px;
                ">Please name the new dataset:</a>
                <el-input 
                    v-model="datasetName"
                    placeholder="Dataset Name"
                    style="
                        width:350px;
                        margin-bottom:30px;
                        margin-top:10px;
                    ">
                <template slot="append">.h5ad</template>
                </el-input>
                <div
                    style="
                        width:350px;
                        display:flex;
                        justify-content:space-between
                    ">
                    <el-button
                        type="success"
                        style="width:150px"
                        @click="Confirm">Confirm</el-button>
                    <el-button
                        type="danger"
                        style="width:150px"
                        @click="Cancel">Cancel</el-button>
                </div>
            </div>
        </el-dialog>
</template>

<script>
import axios from "axios";
import { saveAs } from 'file-saver';

import {createNewJob,loadExistJob} from "@/utils/interface"
import { Form, Tabs, TabPane, Loading, Drawer,Message } from "element-ui";
import {saveLocalDataset,clearCallback} from '@/utils/interface'

export default {
    name: "LocalDatasetPanel",

    props:[

    ],

    watch:{

    },

    data(){
        return {
            datasetName:'',
            showDialog:false,
        }
    },

    computed:{
        curData(){
            return this.$store.state.curData;
        },
        JobId(){
            return this.$store.state.JobId;
        },
        ViewId(){
            return this.curData.ViewId
        },
        chosenData(){
            return this.$store.state.curData.chosenData;
        },
    },

    methods:{
        /**
         * 
         * 内外api
         * 
         */
        openDialog(){//显示对话框
            this.showDialog = true;
        },
        closeDialog(){//关闭对话框
            this.showDialog = false;
        },


        /**
         * 
         * 内部api
         * 
         */
        Confirm(){
            //条件检验1：不允许空名称
            if(this.datasetName == '' || this.datasetName === undefined || this.datasetName === null){
                this.$message({
                    message:'Can not use empty dataset name',
                    type:'error',
                    'showClose':true,
                })
                return;
            }
            //条件检验2：名称不允许包含非法字符
            let forbiddenCharas = ['/','+','-','.']
            for(let c of forbiddenCharas){
                if(this.datasetName.includes(c)){
                    this.$message({
                        'message':`Invalid name: contains forbidden character "${c}"`,
                        'type':'error',
                        'showClose':true,
                    })
                    return;
                }
            }


            const loading = Loading.service({ fullscreen: true });
            loading.text = `Exporting, please wait...`
            saveLocalDataset(this.JobId,this.ViewId,this.chosenData)
                .then((response)=>{
                    const file = response.data
                    let blob = new Blob([file]);
                    // let a = document.createElement('a');
                    // let url = window.URL.createObjectURL(blob);

                    // a.href = url;
                    // a.download = `${this.datasetName}.h5ad`;
                    // a.click();
                    // window.URL.revokeObjectURL(url);
                    saveAs(blob, `${this.datasetName}.h5ad`);

                    
                    this.closeDialog();
                    this.datasetName = '';

                    loading.close();

                    clearCallback(this.JobId,'api/saveLocalDataset',{'ViewId':this.ViewId}).then(res=>console.log(res)).catch(err=>console.log(err))


                })
                .catch((err)=>{
                    console.log(err)
                    this.$message({
                        message:'Save action failed',
                        type:'error',
                        'showClose':true,
                    })

                    this.closeDialog();
                    this.datasetName = '';

                    loading.close()

                    clearCallback(this.JobId,'api/saveLocalDataset',{'ViewId':this.ViewId}).then(res=>console.log(res)).catch(err=>console.log(err))

                })

        },
        Cancel(){
            this.closeDialog();
        },
    },
}
</script>

<style scoped lang="less">
.localDatasetContainer{
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 500px;
        height: 250px;
        // background-color: #24292f;
    }
/deep/ .el-dialog__body{//保存子数据集对话框的主体（除去标题）
    background-color:#F5F5F5;
    border-top: 2px solid lightgray;
    padding: 10px 10px;
}
</style>