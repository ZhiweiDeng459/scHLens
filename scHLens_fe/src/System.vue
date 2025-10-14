<template>
    <div class="content">
        <!--job加载框-->
        <el-dialog
            :visible.sync="unfinishLogIn"
            width="590px"
            :show-close="false"
            :close-on-click-modal="false"
            :close-on-press-escape="false"
            :lock-scroll="true"
            :modal="false"
            >
            <LogInPanel 
                :closeLogIn="closeLogIn"></LogInPanel>
        </el-dialog>
        <!--加载job前遮罩-->
        <div
            class="dialog-modal"
            v-if="unfinishLogIn">
        </div>
        <!--服务器未启动连接提示框-->
        <div 
            style="
                position: absolute;
                z-index:9995;
                height:100px;
                background-color:white;
                top: 40%;
                left: 50%;
                transform: translate(-50%, -50%);
                opacity: 1;
                border-radius:4px;
                box-shadow:0 1px 3px rgba(0,0,0,.3);
                border:4px solid black"
            v-if="!server_connect_info_show">
            <div style="width:100%;height:100%;display: flex;flex-direction: column;align-items: center;justify-content: space-around;">
                <b style="font-size:22px;font-family: YaHei;margin: 0px 20px;">{{ server_connect_text }}</b>
                <b>please wait...</b>
            </div>
        </div>
        <!--服务器未启动遮罩-->
        <div
            class="server-start-modal"
            v-if="!server_connect_info_show">
        </div>

        <!--视图布局-->
        <div class="left">
            <ControlPanel class="control" :layer="curData.TreeNode.layer"></ControlPanel>
        </div>
        <div v-show="dataList.length != 0" class="right">  
            <div class="top">
                    <div class="c-info">
                        <div class="baseinfo-shower-container">
                            <el-input v-model="curCellNum" :readonly="true" size="mini">
                                <template slot="prepend">
                                    <a style="color:white;font-size:15px">Cells</a>
                                </template>
                            </el-input>
                            <el-input v-model="curGeneNum" :readonly="true" size="mini" style="margin-left: 5px;">
                                <template slot="prepend">
                                    <a style="color:white;font-size:15px">Genes</a>
                                </template>
                            </el-input>
                        </div>
                        <div class="gene-shower-container">
                            <div style="background-color: #409EFF;border:1px solid #409EFF;border-radius: 4px;border-top-right-radius:0;border-bottom-right-radius:0;padding: 4px 20px;">
                                <a style="color:white;font-size:15px">Current Gene</a>
                            </div>
                            <el-select v-model="curGeneName" multiple style="width:230px" size="mini" @remove-tag="handleRemoveCurGeneFromSelect" :collapse-tags="true">
                                <el-option
                                    v-for="gene in curGeneName"
                                    :key="gene"
                                    :label="gene"
                                    :value="gene"
                                    style="padding: 0px;"> 
                                    <div style="width: 100%;height: 100%;background-color: white;display: flex;align-items: center;justify-content: space-between;">
                                        <span style="float: left;margin-left: 13px;">{{ gene }}</span>
                                        <i class="el-icon-close" style="color:black;float: right;margin-right: 13px;" @click="handleRemoveCurGeneFromSelect2(gene)"></i>
                                    </div>

                                </el-option>
                            </el-select>
                            <el-button style="border-top-left-radius:0;border-bottom-left-radius:0;" type="danger" size="mini" @click="resetCurGenes">Reset</el-button>
                        </div>
                        <div style="flex:1 1 0">
                            <!--填充块-->
                        </div>
                        <div class="gene-chooser-container">
                            <el-autocomplete
                                v-model="tempCurGeneName"
                                size="mini"
                                :fetch-suggestions="geneQuerySearch"
                                @select="handleGeneSelect"
                                style="width:300px;"
                                placeholder='Search gene to add to "Current Gene"'
                                >
                                <template  slot="prepend">
                                    <!-- <a style="color:white;font-size:15px">Current Gene</a> -->
                                    <i class="el-icon-search" style="color:white"></i>
                                    </template>
                                <!-- <div slot="append">
                                    <el-popover
                                        placement="right"
                                        trigger="click"
                                        class="recommendGenePopoverInGP">
                                            <el-table :data="recommendGenes"  @row-click="selectRecommendGene" style="width:100%" height="400px">
                                                <el-table-column class-name="recommend-item" width="200" property="gene" label="Recommended Genes" align="center"></el-table-column>
                                            </el-table>
                                            <el-button  slot="reference" type="warning" @click="getRecommendGene" icon="el-icon-star-on"></el-button>
                                    </el-popover>
                                </div> -->
                            </el-autocomplete> 
                            <!-- <el-button style="border-top-left-radius:0;border-bottom-left-radius:0;" slot="reference" type="primary" size="mini" @click="handleAddGeneClick" icon="el-icon-right"></el-button> -->
                            <el-popover
                                placement="bottom-start"
                                width="300"
                                title="Multiple genes input"
                                v-model="showMultiGenesInput"
                                trigger="click">
                                <div style="display: flex;flex-direction: column;align-items: center;">
                                    <el-input
                                        type="textarea"
                                        :autosize="{minRows:4,maxRows:4}"
                                        v-model="multiGenesText"
                                        placeholder="You can input multiple genes separated by spaces, line breaks, or semicolons here. Then scHLens will push them into the 'current gene'."></el-input>
                                    <div
                                        style="
                                            width:280px;
                                            display:flex;
                                            justify-content:space-between;
                                            margin-top:15px;
                                        ">
                                        <el-button
                                            type="primary"
                                            style="width:120px"
                                            size="mini"
                                            @click="onMultiGenesConfirm">Confirm</el-button>
                                        <el-button
                                            type="danger"
                                            style="width:120px"
                                            size="mini"
                                            @click="onMultiGenesCancel">Cancel</el-button>
                                    </div>            
                                </div>
                                <el-button slot="reference" style="border-top-left-radius:0;border-bottom-left-radius:0;" type="primary" size="mini" icon="el-icon-s-unfold"></el-button>
                            </el-popover>
                        </div>
                        <el-button @click="openGroupPanel" size="mini" type="primary">Group</el-button>
                        <el-button @click="refreshAll" size="mini" type="primary" icon="el-icon-refresh"></el-button>
                    </div>
                    <div class="top-views-container">
                        <view-container style="flex:0 0 500px;margin-right:5px" :type="'CellProjection'"/>
                        <view-container style="flex:0 0 500px;margin-right:5px" :type="'GeneProjection'"/>
                        <view-container style="flex:1 1 0" :type="'GeneExpression'"/>
                    </div>
            </div>
            <div class="bottom">
                <view-container style="flex:0 0 600px;margin-right: 5px;" :type="'CellChat'"/>
                <view-container style="flex:1 1 0;margin-right:5px" :type="'MarkerGene'"/>
                <ProjectionTree style="flex:0 0 400px;min-width:0px;"/>
            </div>
        </div>

        <!--GroupPanel-->
        <GroupPanel ref="group-panel"></GroupPanel>
        <!--无数据遮罩层-->
        <div v-show="dataList.length == 0" class="empty-modal-container">
            <el-empty description="NO DATA" :image-size="500" style="height:100%;width:100%">
            </el-empty>
        </div>
        

    </div>
</template>

<script>
import Vue from "vue";
import {Step,Steps,Message,Empty,Drawer,Notification} from "element-ui";

import {InstanceClose,recommendGene,requestCandidateGeneList,multiGenesSplitFromText,testConnection} from '@/utils/interface.js';
import ControlPanel from "@/views/ControlPanel";
import ViewContainer from "@/components/ViewContainer"
// import ViewChooser from "@/components/ViewChooser"
import CellProjection from "@/views/CellProjection";
import GeneProjection from "@/views/GeneProjection";
import GeneExpression from "@/views/GeneExpression";
import MarkerGene from "@/views/MarkerGene";
import ProjectionTree from "@/views/ProjectionTree";
import TrajectoryInference from  "@/views/TrajectoryInference"
import LogInPanel from "@//components/LogInPanel"
import GroupPanel from "@/components/GroupPanel"
import MessageBoard from "./components/MessageBoard.vue";

Vue.component(Empty.name,Empty)
Vue.component(Step.name, Step);
Vue.component(Steps.name, Steps);
Vue.component(Drawer.name,Drawer)
Vue.component(Notification.name,Notification)

Vue.prototype.$message = Message;
export default {
    name: "System",
    components: {
        ControlPanel,
        ViewContainer,
        LogInPanel,
        ProjectionTree,
        GroupPanel,
    },
    computed: {
        curData() {
            return this.$store.state.curData;
        },
        JobId(){
            return this.$store.state.JobId;
        },
        ViewId(){
            return this.curData.ViewId
        },
        curGeneName(){
            return this.$store.state.curGeneName;
        },
        genes(){
            return this.curData.genes
        },
        dataList(){
            return this.$store.state.dataList
        },
        curCellNum(){
            if(this.$store.state.curData === undefined || this.$store.state.curData === null){
                return 0
            }
            if(this.$store.state.curData.cellData === undefined || this.$store.state.curData.cellData === null){
                return 0
            }
            else{
                return this.$store.state.curData.cellData.length;
            }
        },
        curGeneNum(){
            if(this.$store.state.curData === undefined || this.$store.state.curData === null){
                return 0
            }
            if(this.$store.state.curData.genes === undefined || this.$store.state.curData.genes === null){
                return 0
            }
            else{
                return this.$store.state.curData.genes.length;
            }

        },
        messageBoard(){
            return this.$store.state.messageBoard
        }
    },
    methods:{
        closeHandler(event){ //网页关闭或者刷新
            InstanceClose()
        },
        closeLogIn(){//关闭登录窗
            this.unfinishLogIn = false;
        },
        geneQuerySearch(queryString, cb){//基因选择器的候选搜索


            // requestCandidateGeneList(this.JobId,this.curData.ViewId,"^" + queryString)
            //     .then((response) => {
            //         const data = response.data;
            //         cb(
            //             data.map((item) => {
            //                 return { value: item };
            //             }).reverse()
            //         );
            //     })
            //     .catch((err) => {
            //         console.log(err);
            //     });

            let matched_genes = this.curData.queryGenes.filter(str => str.toLowerCase().startsWith(queryString.toLowerCase()))
            let sorted_matched_genes = matched_genes.sort((a,b)=>a.length - b.length)
            cb(sorted_matched_genes.map(item=>{
                return { value: item };
            }))
            
        },
        handleGeneSelect(item){//处理基因选择（新增基因进入curGeneName）
            if(this.curGeneName.includes(item.value)){//基因已经有了该基因
                this.$message({
                    'message':'This Gene has been added',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            this.$store.commit("addToCurGeneName", item.value);

            
        },
        // getRecommendGene(){ //获取推荐基因
        //     this.recommendGenes.length = 0;
        //     recommendGene(this.JobId,this.curData.ViewId,this.$store.state.recommendMode)
        //     .then((response) => {
        //         for(let item of response.data){
        //             this.recommendGenes.push({
        //                 'gene':item
        //             })
        //         }
        //     })
        //     .catch((err) => {
        //         console.log(err);
        //     });
        // },
        // selectRecommendGene(row){//点击了推荐基因表格中的基因（新增基因进入curGeneName）
        //     if(this.curGeneName.includes(row.gene)){//基因已经有了该基因
        //         this.$message({
        //             'message':'This Gene has been added',
        //             'type':'error',
        //             'showClose':true,
        //         })
        //         return;
        //     }
        //     this.$store.commit("addToCurGeneName", row.gene);
        // },
        reDrawTree(){
            this.$store.$emit("refreshView","ProjectionTree")
        },
        refreshAll(){//刷新所有视图
            this.$store.commit("refreshView",'all')
            // this.$message({
            //     'message':'All Views has been updated',
            //     'type':'success',
            //     'showClose':true,
            // })
        },
        openGroupPanel(){//打开group视图
            this.$refs['group-panel'].openDialog();
            this.r_Drawer_visible = true
        },
        handleRemoveCurGeneFromSelect(value){//对于当前值显示框中的移除按钮
            if(this.curGeneName.length<=1){
                this.$message({
                    'message':'Retain at least one gene',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            this.$store.commit("deleteFromCurGeneName", value);
        },
        handleRemoveCurGeneFromSelect2(value){//对于下来框中的删除按钮
            if(this.curGeneName.length<=1){
                this.$message({
                    'message':'Retain at least one gene',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            this.$store.commit("deleteFromCurGeneName", value);
        },
        // handleAddGeneClick(){//点击了基因搜索框中确定添加基因的按钮

        // },
        // openMultiGenesInput(){//打开多基因输入框
        //     this.showMultiGenesInput = true;
        // },
        resetCurGenes(){//重设currentGenes
            this.$store.commit("initCurGeneName");
        },
        onMultiGenesConfirm(){//多基因输入确认
            multiGenesSplitFromText(this.JobId,this.ViewId,this.multiGenesText)
                .then((response)=>{
                    let data = response.data;
                    let split_num = data.split_num;
                    let valid_num = data.valid_num;
                    let invalid_num = data.invalid_num;
                    let valid_genes = data['valid_genes'];//确保每个基因都是合法的，即在数据中已经出现的
                    let invalid_genes = data['invalid_gens']

                    //提示
                    this.messageBoard.show()
                    let message = {
                        'Overview':`Among the ${split_num} detected genes, ${valid_num} are valid(included in the dataset), while ${invalid_num} are invalid(not included in the dataset).`,
                        'Invalid Genes':invalid_num == 0 ? '-' :invalid_genes.join(" "),
                        'Valid Genes':valid_num == 0 ? '-' : valid_genes.join(" "),
                    }
                    this.messageBoard.setTitle('Multiple genes input')
                    this.messageBoard.setMessageData(message)

                    //上传
                    let union_genes = [...new Set([...valid_genes, ...this.curGeneName])]
                    this.$store.commit("updateCurGeneName", union_genes);
                    


                }).catch((err)=>{
                    this.$message({
                        'message':'Fail to split',
                        'type':'error',
                        'showClose':true,
                     })
                })
                this.multiGenesText = '';
                this.showMultiGenesInput = false;
        },
        onMultiGenesCancel(){//多基因输入取消
            this.multiGenesText = '';
            this.showMultiGenesInput = false;
        }
    },

    watch:{
        curGeneName:{
            deep:true,
            handler(){
                //提示基因已经切换
                // this.$message({
                //     'message':'Current Gene has been updated',
                //     'type':'success',
                //     'showClose':true,
                // })
            }
        }

    },

    data(){
        return{
            unfinishLogIn:true,
            tempCurGeneName:'',
            //recommendGenes:[],
            r_Drawer_visible:false,
            showMultiGenesInput:false,//展示多基因输入框
            multiGenesText:'',//输入的多基因文本

            server_connect_info_show:false, //是否显示后台服务连接的提示
            server_connect_text : 'The backend server is starting...', //后台服务连接显示的文本
        }
    },
    created(){
        window.addEventListener('beforeunload',this.closeHandler)

        //连接检测：
        testConnection().on('test',data=>{
            if(data['status'] == 0){
                //后端已经启动，关闭遮罩层
                this.server_connect_info_show = true
            }
        }).on('disconnect',data=>{//断开连接后恢复加载页面

            // this.server_connect_info_show = false;
            // this.server_connect_text = 'Connection lost. Reconnecting to the server...';
        })

    },

};
</script>

<style scoped lang="less">
.content {
    flex: 1 1 0;
    display: flex;
    align-items: stretch;
    position: relative;
    background-color:#F5F5F5;
    overflow-x:hidden;
    .left {
        flex: 0 0 220px;
        display: flex;
        flex-direction: column;
        align-items: center;
        margin: 0px 0px 0px 0px;
        background-color:#24292f;
        border-right: 2px solid rgb(200, 200, 200);
        .SystemTitle{
            font-size: 50px;
            font-family: AslinaBold;
            color: #90e36b;
            margin-top:15px;
            margin-bottom:20px;
        }
        .control {
            margin-top:30px;
            flex: 1 1 0;
        }
    }
    .empty-modal-container{
        flex: 1 1;
        position: relative;  
        /deep/ .el-empty__description p{
            width:500px;
            font-size:30px
        }
    }
    .right {
        flex: 1 1;
        position: relative;
        display: flex;
        flex-direction: column;
        align-items: stretch;
        min-width: 0;
        margin: 5px 10px 7px 10px;
        .top{
            display: flex;
            margin-bottom:5px;
            flex: 0 0 500px;
            flex-direction: column;
            .c-info{
                display: flex;
                flex:0 0;
                align-items: center;
                margin: 0px 0px 3px 2px;
                .baseinfo-shower-container{//基础信息演示容器
                    display: flex;
                    align-items: center;
                    margin: 0 7px 0px 0px;
                    /deep/ .el-input-group{
                        
                        .el-input-group__prepend{
                            background-color: #323232 !important;
                            border:1px solid #323232;
                        }
                        // .el-input__inner{
                        //     // border-top:1px solid #666666;
                        //     // border-bottom:1px solid #666666;
                        //     border:1px solid #666666;
                        // }
                    }

                }
                .gene-shower-container{
                    margin: 0 7px 0px 0px;
                    display: flex;
                    align-items: center;

                }
                .gene-chooser-container{//自定义的基因选择器
                    display: flex;
                    align-items: center;
                    background-color: white;
                    margin: 0 7px 0px 0px;
                    // border-top: 2px solid rgb(200, 200, 200);
                    // border-left: 2px solid rgb(200, 200, 200);
                    // border-right: 2px solid rgb(200, 200, 200);
                    // border-radius: 2px;
                    /deep/ .el-input-group{
                        
                        .el-input-group__prepend{
                            background-color: #409EFF !important;
                            border:1px solid #409EFF;
                        }
                        // .el-input__inner{
                        //     // border-top:1px solid #666666;
                        //     // border-bottom:1px solid #666666;
                        //     border:1px solid #666666;
                        // }
                    }

                    /deep/ .el-input-group__append{
                            background-color: #E6A23C !important;
                            border:1px solid #E6A23C;
                            .el-icon-star-on{
                                color:white;
                            }
                    }

                }
            }

            .top-views-container{
                display: flex;
                flex: 1 1 0;
                flex-direction: row;
                align-items: stretch;
            }

        }
        .bottom{
            flex: 1 1 0;
            min-height: 0;
            display: flex;
            flex-direction: row;
            align-items: stretch;
            background-color: white;
        }


    }





}

//对话框-遮罩层
.el-dialog__wrapper{
    position: absolute;
}
.server-start-modal{
    position: absolute;
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    z-index: 9994;
    opacity: 0.6;
    background: #000;
  
}

.dialog-modal{
    position: absolute;
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    z-index: 1;
    opacity: 0.5;
    background: #000;
}





/deep/ .recommend-item{
    cursor: pointer;
}

.el-input.is-readonly .el-input__inner {
    color: black; /* 设置字体颜色为灰色 */
}

.el-textarea__inner::placeholder {
  white-space: pre-line;      /* 保留换行符，但允许换行 */
  word-break: keep-all;       /* 不要把单词拆开 */
  overflow-wrap: break-word;  /* 单词太长时才强制断行 */
}

</style>