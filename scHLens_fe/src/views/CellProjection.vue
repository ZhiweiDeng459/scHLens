<template>
        <div class="cell-projection-container">
            <div class="view-header">
                <div style="flex:1 1 0">
                    <a id="cell-projection-title">Cell Projection View</a>
                </div>
            </div>
            <scatter ref="scatter" class="scatter-view" v-show="activeFlag"></scatter>
        </div>
</template>

<script>
import Scatter from "@/components/CellProjection/Scatter";
import { Loading } from "element-ui";
import * as d3 from "d3";
import { generateViewID } from "@/utils/IDgenerator";
import eventBus from "@/utils/eventBus.js"


export default {
    name: "CellProjection",
    components: {
        Scatter,
    },
    computed: {
        curData() {
            return this.$store.state.curData;
        },
        cellData() {
            return this.curData.cellData;
        },
        chosenData() {
            return this.curData.chosenData;
        },
        activeFlag() {
            return this.curData.activeFlag['CellProjection'];
        },
        JobId(){
            return this.$store.state.JobId;
        }
    },

    methods: {
    },

    mounted(){
        //eventbus控制loading
        let loadingArray = [];
        eventBus.$on('CellProjectionViewRefreshingStart',()=>{
            loadingArray.push(Loading.service({
                target:".scatter-view",
                lock:true,
                text:"Refreshing",
                spinner: 'el-icon-loading',
                background: 'rgba(255, 255, 255, 0.8)',
            }));
        })
        eventBus.$on('CellProjectionViewRefreshingClose',()=>{
            for(let loading of loadingArray)
                loading.close();
        })
    }

};
</script>

<style scoped lang="less">
.cell-projection-container {
    position: relative;
    width: 100%;
    height: 100%;
    display: flex;
    flex-direction: column;
    align-items: stretch;
    // border: 2px solid rgb(200, 200, 200);
    // border-radius: 10px;
    .view-header{
        padding-left:5px;
        flex:0 0 40px;
        display: flex;
        align-items: center;
        background-color: #24292f;
        border:2px solid #24292f;
        #cell-projection-title{
            // font-family:YaHei;
            color:white;
            font-size:20px;
        }   
    }
    .scatter-view {
        flex: 1 1 0;
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
