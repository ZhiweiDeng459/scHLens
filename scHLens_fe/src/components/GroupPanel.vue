<template>
        <div>
        <el-dialog
            :visible.sync="showDialog"
            width="650px"
            :show-close="true"
            :lock-scroll="false"
            :modal="false"
            :close-on-click-modal="false"
            title="Group Panel"
            v-dialogDrag
            >
            <div slot="title">
                <b style="font-size:20px;">Group Info</b>
            </div>
            <div class="GroupInfoContainer">
                <div>

                </div>
                <el-scrollbar class="group-scroller" >
                    <div v-for="(group,index) in localGroups" :key="`${JobId}-${ViewId}-${group.id}`" class="group-table-row">

                        <div style="display:flex;align-items:center">
                            <el-color-picker popper-class="hhh" v-model="group.color" @change="handleColorChangeFinish(group.id)"></el-color-picker>
                            <el-tooltip class="item" effect="dark" :content="group.id" placement="top">
                                <b style="margin-left:15px;font-family:YaHei;font-size:16px">{{group.id.length>5?group.id.substr(0,5) + '...':group.id}}</b>
                            </el-tooltip>
                            <a style="margin-left:10px;color:gray;margin-top:2px">{{`(${group.size})`}}</a>
                        </div>



                        <div>
                            <el-tooltip class="item" effect="dark" :content="'Export Marker Genes'" placement="top">
                                <el-button icon="el-icon-download" @click="handleExportMarkerGenes(group.id)" size="mini" style="width:26px;height:26px;padding:0px;margin-right:10px;" type="primary"></el-button>
                            </el-tooltip>

                            <el-popover
                                placement="right"
                                width="700"
                                @show="handleAnnoRecomShow(index)">
                                <el-input @keyup.enter.native="handleAnnotationFinish(group.id)" @blur="handleAnnotationFinish(group.id)" slot="reference" v-model="group.name" size="mini" style="width:350px">
                                    <template slot="prepend">Annotation</template>
                                </el-input>
                                <AnnoRecom @curGeneSetChange="handleAnnoRecomCurGeneSetChange" @RecommendationChosen="handleRecommendationChosen" :ref="`AnnoRecom`" :cluster_id="group.id"/>
                            </el-popover>
                        </div>
                    </div>
                </el-scrollbar>
            </div>
        </el-dialog>
        </div>
</template>

<script>
import Vue from 'vue'
import axios from "axios";
import {createNewJob,loadExistJob,updateGroupName,updateGroupColor,exportGlobalMarkers,exportLocalMarkers} from "@/utils/interface"
import { Form, Tabs, TabPane, Loading, Drawer,Message,ColorPicker} from "element-ui";
import {saveLocalDataset,clearCallback} from '@/utils/interface'
import AnnoRecom from '@/components/AnnoRecom'

Vue.component(ColorPicker.name,ColorPicker);

export default {
    name: "GroupPanel",
    components:{
        AnnoRecom
    },
    props:[

    ], 

    data(){
        return {
            showDialog:false,
            localGroups:[],//groups的本地内容，后续更新
            AnnoRecomVisible:[],

            //AutoRecommen
            curGeneSet:[],//当前的物种类型和基因集 ['Human', 'PanglaoDB Augmented 2021']
            
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
        groups(){
            return this.curData.groups
        }
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
        commitLocalGroupName(cluster_id){//将当前localGroup(cluster_id指定)的name更新到全局groups和后端
            let group = this.groups.find(g=>g.id == cluster_id)
            let localGroup = this.localGroups.find(g=>g.id == cluster_id)
            group.name = localGroup.name
            //group_name修改结果提交给后端
            let new_group_name = {}
            for(let g of this.groups){
                new_group_name[g.id] = g.name
            }
            updateGroupName(this.JobId,this.ViewId,new_group_name)
                .then((response)=>{
                })
                .catch((err) => {
                    console.log(err);
                });
        },
        commitLocalGroupColor(cluster_id){//将当前localGroup(cluster_id指定)更新到全局groups和后端
            let group = this.groups.find(g=>g.id == cluster_id)
            let localGroup = this.localGroups.find(g=>g.id == cluster_id)
            group.color = localGroup.color
            console.log('color:',localGroup.color)
            //group_color修改结果提交给后端
            let new_group_color = {}
            for(let g of this.groups){
                new_group_color[g.id] = g.color
            }
            updateGroupColor(this.JobId,this.ViewId,new_group_color)
                .then((response)=>{

                })
                .catch((err)=>{
                    console.log(err)
                })

        },
        


        /**
         * 
         * handle
         * 
         */
        handleAnnoRecomShow(index){//推荐框弹出，弹出时将最近的一次基因集选择记录传入
            this.$refs['AnnoRecom'][index].show(this.curGeneSet);//同时会将所有推荐窗口中，最新选择的基因库传递给弹出的推荐窗口
        },
        handleAnnotationFinish(cluster_id){ //用户手动完成了一次注释，提交修改后的注释到groups，并向后端更新
            this.commitLocalGroupName(cluster_id);
        },  
        handleColorChangeFinish(cluster_id){ //用户修改了一次聚类颜色，提交修改后的颜色给groups，并向后端更新
            //检查修改颜色的合法性
            let group = this.groups.find(g=>g.id == cluster_id)
            let localGroup = this.localGroups.find(g=>g.id == cluster_id)
            if(localGroup.color===null){//如果点击了clear
                localGroup.color = group.color;
                return;
            }
            else if(localGroup.color==''){//如何输入了不合法的颜色
                localGroup.color = group.color;
                this.$message({
                    'message':'Please input legal color such as "#FFFFFF"',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            //提交
            this.commitLocalGroupColor(cluster_id)
        },
        handleRecommendationChosen(params){//用户从标注框内选择了一个细胞类型，提交细胞类型给对应的localgroup，更新到全局group和后端中
            let localGroup = this.localGroups.find(g=>g.id == params.cluster_id)
            localGroup.name = params.type;
            this.commitLocalGroupName(params.cluster_id)
        },
        handleAnnoRecomCurGeneSetChange(curGeneSet){//用户在推荐框切换了新的gene set，这里保存用于最近使用的gene set的记录
            this.curGeneSet = curGeneSet;
        },
        handleExportMarkerGenes(groupId){//导出某一个聚类的marker genes
            //分别导出全局基因和局部基因
            const self = this;
            //首先检查是全局还是局部
            let layer = self.curData.TreeNode.layer
            //根据groupId检索groupname
            let group_name = self.groups.find(g=>g.id==groupId)['name']

            const loading = Loading.service({ fullscreen: true });

            if(layer != 1){//如果不为全局节点，那么导出局部marker
                exportLocalMarkers(self.JobId,self.curData.ViewId,groupId)
                    .then((response)=>{
                        const file = response.data
                        let blob = new Blob([file]);
                        let a = document.createElement('a');
                        let url = window.URL.createObjectURL(blob);
                        a.href = url;
                        a.download = `LocalMarkers_${self.JobId}_${self.curData.ViewId}_${group_name}.xlsx`;
                        a.click();
                        window.URL.revokeObjectURL(url);
                        clearCallback(this.JobId,'api/exportLocalMarkers',{'ViewId':this.ViewId}).then(res=>console.log(res)).catch(err=>console.log(err))

                    })
                    .catch((err)=>{
                        console.log(err)
                        this.$message({
                            'message':'An error occurred during export. This error may be due to not selecting the marker procedure in the pipeline, or other unknown reasons.',
                            'type':'error',
                            'showClose':true,
                        })
                        clearCallback(this.JobId,'api/exportLocalMarkers',{'ViewId':this.ViewId}).then(res=>console.log(res)).catch(err=>console.log(err))

                    })
            }
            //导出全局marker
            exportGlobalMarkers(self.JobId,self.curData.ViewId,groupId)
                    .then((response)=>{
                        const file = response.data
                        let blob = new Blob([file]);
                        let a = document.createElement('a');
                        let url = window.URL.createObjectURL(blob);
                        a.href = url;
                        a.download = `GlobalMarkers_${self.JobId}_${self.curData.ViewId}_${group_name}.xlsx`;
                        a.click();
                        window.URL.revokeObjectURL(url);
                        loading.close()
                        clearCallback(this.JobId,'api/exportGlobalMarkers',{'ViewId':this.ViewId}).then(res=>console.log(res)).catch(err=>console.log(err))

                    })
                    .catch((err)=>{
                        console.log(err)
                        this.$message({
                            'message':'An error occurred during export. This error may be due to not selecting the marker procedure in the pipeline, or other unknown reasons.',
                            'type':'error',
                            'showClose':true,
                        })
                        loading.close()
                        clearCallback(this.JobId,'api/exportGlobalMarkers',{'ViewId':this.ViewId}).then(res=>console.log(res)).catch(err=>console.log(err))

                    })

        }

    },

    watch:{
        groups:{//时刻用group的最新去更新localgroup
            deep:true,
            handler(){
                if(this.groups === undefined || this.groups === null)
                    this.localGroups = []
                else{
                    this.localGroups = JSON.parse(JSON.stringify(this.groups))
                }
            }

        }
    },
    

}
</script>

<style scoped lang="less">
.GroupInfoContainer{//容器
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 100%;
        height: 100%;
        // background-color: #24292f;

}

.group-header{//表头

    width: 100%;
    display: flex;

}

.group-scroller{//滚动层
    width: 100%;
    height: 100%;
}


.group-table-row{//表行
    margin:5px 10px;
    padding: 10px 10px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    border: 2.5px solid #e5e5e5;
    border-radius: 7px;
}




/**
* scroll
*/

/deep/ .el-scrollbar__wrap{
    overflow: hidden;
    position: relative;
    width:100%;
}
/deep/ .is-horizontal{
    height: 10px;
    .el-scrollbar__thumb{
        background-color:rgb(150, 150, 150);
    }
}
/deep/ .is-vertical{
    width: 10px;
    .el-scrollbar__thumb{
        
        background-color:rgb(150, 150, 150);
    }
}


/**
* dialog
*/
/deep/ .el-dialog{
    pointer-events: auto;
    border: 2px solid lightgray;
    
    .el-dialog__body{//对话框的主体（除去标题）
        background-color:white;
        padding: 10px 10px;
        height: 300px;//决定了dialog的高度
    }
    .el-dialog__header{//对话框的标题
        background-color:rgba(0, 0, 0, 0.03);
        border-bottom: 1px solid rgba(0, 0, 0, 0.125);
        .el-dialog__title{
            color:black;
            font-size:20px;
        }
        .el-dialog__close{
            color:black
        }
    }

}
.el-dialog__wrapper {
    pointer-events: none;
  }

</style>