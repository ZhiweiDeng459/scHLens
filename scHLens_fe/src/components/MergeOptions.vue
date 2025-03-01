<template>
        <el-dialog
            :visible.sync="showDialog"
            width="525px"
            :show-close="true"
            :lock-scroll="true"
            :modal="true"
            :close-on-click-modal="false"
            top="30vh"
            >
            <div slot="title">
                <b style="font-size:20px;">Merge Views Settings</b>
            </div>
            <div class="MergeOptionsContainer">
                <div class="option-row">
                    <b class="option-title">
                        Merge with Projection
                        <el-tooltip style="margin-left:2px;" content='Merge the projection effect of the local node onto the global node.' placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </b>
                    <div class="option-content">
                        <el-switch
                            v-model="mergeOption['isProjection']"></el-switch>
                    </div>
                </div>
                <div class="option-row">
                    <b class="option-title">
                        Merge with Labels(Annotations)
                        <el-tooltip style="margin-left:2px;" content='Replace the corresponding labels in the global node with the labels from the local node when merging views.' placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </b>
                    <div class="option-content">
                        <el-switch
                            v-model="mergeOption['isLabel']"></el-switch>
                    </div>
                </div>
                <div class="option-row">
                    <b class="option-title">
                        Merge Labels(Annotations) with the same name
                        <el-tooltip style="margin-left:2px;" content="Merge labels with the same name when the option 'Merge with Labels' is selected" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </b>
                    <div class="option-content">
                        <el-switch
                            v-model="mergeOption['isMergeSameNameLabels']"></el-switch>
                    </div>
                </div>
                <div class="option-row">
                    <b class="option-title">
                        Merge small Labels(Annotations)
                        <el-tooltip style="margin-left:2px;" content="Cells with labels that include too few cells will be merged with the labels of the nearest clusters; otherwise, it might not be possible to calculate marker genes." placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </b>
                    <div class="option-content">
                        <el-switch
                            v-model="mergeOption['isMergeSmallLabels']"></el-switch>
                    </div>
                </div>


                <div
                    style="
                        width:400px;
                        display:flex;
                        justify-content:space-between;
                        margin-top:10px;
                    ">
                    <el-button
                        type="primary"
                        size="small"
                        style="width:150px"
                        @click="Confirm">Confirm</el-button>
                    <el-button
                        type="danger"
                        size="small"
                        style="width:150px"
                        @click="Cancel">Cancel</el-button>
                </div>
            </div>
        </el-dialog>
</template>

<script>
import Vue from "vue"
import axios from "axios";
import {createNewJob,loadExistJob} from "@/utils/interface"
import {CheckboxGroup,Checkbox,Form,FormItem,Switch} from "element-ui";
import {saveLocalDataset} from '@/utils/interface'
Vue.component(CheckboxGroup.name, CheckboxGroup);
Vue.component(Checkbox.name, Checkbox);
Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Switch.name, Switch);

export default {
    name: "MergeOptions",

    props:[
        'mergeFunction'
    ],

    watch:{

    },

    data(){
        return {
            datasetName:'',
            showDialog:false,
            mergeOption:{
                'isLabel':true,
                'isProjection':false,
                'isMergeSmallLabels':true,
                'isMergeSameNameLabels':true,
            }
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
         * 外部api
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
            //check1：merge必须至少选择一项
            if(!this.mergeOption['isLabel'] && !this.mergeOption['isProjection']){
                this.$message({
                    'message':"Please select at least one of 'Merge with Labels' and 'Merge with Labels'.",
                    'type':'error',
                    'showClose':true,
                    })
                return;
            }

            this.mergeFunction(this.mergeOption);
            this.closeDialog();
        },
        Cancel(){
            this.closeDialog();
        },
    },
}
</script>

<style scoped lang="less">
.MergeOptionsContainer{
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 500px;
        height: 250px;
        // background-color: #24292f;
    }
/deep/ .el-dialog__body{//保存子数据集对话框的主体（除去标题）
    background-color:white;
    border-top: 2px solid lightgray;
    padding: 10px 10px;
    height: 240px;//决定了dialog的高度

}

.option-row{
    display: flex;
    flex-direction: row;
    border-bottom: 2px solid lightgray;
    margin-bottom:5px;
    padding:10px 5px;
    width:calc(100% - 30px);
    .option-title{
        flex:0 0 370px;
        font-size:15px;
        font-family: YaHei;  
        color:black;        
    }
    .option-content{
        flex:1 1 0;
        display: flex;
        flex-direction: row;
        justify-content: flex-end;
    }
}

</style>