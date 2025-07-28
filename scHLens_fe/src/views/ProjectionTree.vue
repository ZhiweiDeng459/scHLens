<template>
        <div class="projection-tree-container">
            <div class="view-header">
                <div>
                    <a id="projection-tree-title">Hierarchy Navigation</a>
                </div>
            </div>
            
            <navigation class="navigation-view" :callMergeOptions="callMergeOptions"></navigation>
            <MergeOptions ref="MergeOptions" :mergeFunction="mergeViews"/>
        </div>
</template>

<script>
import Navigation from "@/components/Navigation";
import MergeOptions from "@/components/MergeOptions";
import { mergeViews,fetchViewData } from "@/utils/interface";
import { Loading } from "element-ui";
import eventBus from "@/utils/eventBus.js"


export default {
    name: "ProjectionTree",
    components: {
        Navigation,
        MergeOptions
    },
    computed: {
        dataList() {
            return this.$store.state.dataList;
        },
        projectTree() {
            return this.$store.state.projectTree;
        },
        repaintTag() {
            return this.$store.state.repaintTag;
        },
        JobId(){
            return this.$store.state.JobId;
        },
        mergeRoot(){
            return this.$store.state.mergeRoot
        },
        chosenViews(){
            return this.$store.state.chosenViews;
        },
    },
    methods: {
        mergeViews(mergeOptions) {
            if(this.mergeRoot === undefined){
                this.$message({
                    'message':'Inappropriate set of views was selected!',
                    'type':'error',
                    'duration':10000,
                    'showClose':true,
                })
                return;
            }

            const loading = Loading.service({ fullscreen: true , text:'It may take several minutes, please wait...'});
            let globalViewId = this.mergeRoot
            
            //TODO 这里可能需要注意以下优先级的问题
            let localViewIdList = [];
            for (let view of this.chosenViews) {
                if(view != this.mergeRoot)
                    localViewIdList.push(view);
            }
            mergeViews(this.JobId,globalViewId,localViewIdList,mergeOptions)
                .then((response) => {
                    
                    console.log(response.data)

                    if(response.data.status == 0){//error
                        this.$message({
                            'message':'A error ocurred in running process of the merge process. Maybe you should check whether the local data are proper',
                            'type':'error',
                            'duration':10000,
                            'showClose':true,
                        }) 
                        loading.close();

                    } 
                    else if(response.data.status == 1){//success

                        fetchViewData(this.JobId,this.mergeRoot)
                            .then((fetch_response)=>{
                                let merge_data = fetch_response.data
                                console.log('merge_data',merge_data)
                                this.$message({
                                    'message':'The merge process finished',
                                    'type':'success',
                                    'showClose':true,
                                })
                                this.$store.commit("updateViewData",merge_data)
                                //切换到被合并的视图
                                this.$store.commit("toggleCurData", this.mergeRoot);

                                loading.close();

                            })
                            .catch((fetch_err)=>{
                                console.log(fetch_err)
                                this.$message({
                                    'message':'A error ocurred in running process of the merge process. Maybe you should check whether the local data are proper',
                                    'type':'error',
                                    'showClose':true,
                                })  
                                loading.close();
                            })



                    }
                    else{
                        loading.close();
                    }


                })
                .catch((err) => {
                    console.log(err);
                    loading.close();
                });
        },

        callMergeOptions(){//打开合并视图option面板
            /**
             * 启动选择面板
             */
             this.$refs['MergeOptions'].openDialog();
        }
    },
    mounted(){
        // //eventbus控制loading
        // let loadingArray = [];
        // eventBus.$on('ProjectionTreeRefreshingStart',()=>{
        //     loadingArray.push(Loading.service({
        //         target:".navigation-view",
        //         lock:true,
        //         text:"Refreshing",
        //         spinner: 'el-icon-loading',
        //         background: 'rgba(255, 255, 255, 0.8)',
        //     }));
        // })
        // eventBus.$on('ProjectionTreeRefreshingClose',()=>{
        //     for(let loading of loadingArray)
        //         loading.close();
        // })
    }
};
</script>

<style scoped lang="less">
.projection-tree-container {
    position: relative;
    display: flex;
    flex-direction: column;
    align-items: stretch;
    background-color: white;
    border: 2px solid rgb(200, 200, 200);
    border-radius: 2px;
    .view-header{
        padding:5px;
        flex:0 0 30px;
        display: flex;
        align-items: center;
        background-color: #24292f;
        border:2px solid black;
        #projection-tree-title{
            font-size:20px;
            // font-family:YaHei;
            color:white
        }          
    }
    .el-button-merge {
        text-align: center;
        padding: 5px;
        width: 70px;
        height: 30px;
        .merge-button-font{
            font-size:15px;
        }
    }
    .navigation-view {
        flex: 1 1;
    }
}

//调整loading的字体大小

/deep/ .el-icon-loading{
    font-size:30px;
}

/deep/ .el-loading-text{
    font-size:25px;
}

</style>
