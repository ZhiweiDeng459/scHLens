<template>
        <div class="cell-chat-container">
            <div class="cc-view-header">
                <a id="cell-chat-title">Cell Communication View</a>
                <!-- <el-tooltip content="Display the difference of gene-expression among clusters" placement="top">
                    <i class="el-icon-question"></i>
                </el-tooltip> -->
            </div>
            <div class="cc-view-container">
                <div class="cell-chat-view-container">
                    <CellChatChordPlot v-if="CellChatMode=='Chord'" class="cell-chat-view" :setInteractionTable="setInteractionTable" :dataMode="dataMode"></CellChatChordPlot>
                    <CellChatEdgePlot v-if="CellChatMode=='Edge'"  class="cell-chat-view" :setInteractionTable="setInteractionTable" :dataMode="dataMode"></CellChatEdgePlot>
                </div>
                <div class="interaction-table-container">
                    <InteractionTable style="width:100%;padding-top:40px;" ref="interactionTable"></InteractionTable>
                </div>
                <el-radio-group v-model="CellChatMode" size="mini" class="CellChatMode-toggle">
                    <el-radio-button label="Chord"></el-radio-button>
                    <el-radio-button label="Edge"></el-radio-button>
                </el-radio-group>
                <el-radio-group v-model="dataMode" size="mini" class="dataMode-toggle">
                    <el-radio-button label="weight">Weight</el-radio-button>
                    <el-radio-button label="count">Count</el-radio-button>
                </el-radio-group>
            </div>

        </div>
</template>

<script>
import CellChatChordPlot from "@/components/CellChat/CellChatChordPlot"
import CellChatEdgePlot from "@/components/CellChat/CellChatEdgePlot"
import InteractionTable from "@/components/CellChat/InteractionTable"
import {Loading} from "element-ui";

import eventBus from "@/utils/eventBus.js"

export default {
    name: "CellChat",
    components: {
        CellChatChordPlot,
        CellChatEdgePlot,
        InteractionTable
    },
    data(){
        return {
            CellChatMode:'Chord', //Chord or Edge
            dataMode:'weight', //Weight or Count
        }
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
            return this.curData.activeFlag['CellChat'];
        },
        curMethod(){
            if(this.$store.state.dataList.length == 0){//还没有运行过数据
                return undefined;
            }
            else if(!('CC' in this.$store.state.curData.paramsObj)){//当前节点没有执行过轨迹推断
                return undefined;
            }
            else{
                return Object.keys(this.$store.state.curData.paramsObj.CC)[0]
            }
        }
    },
    methods: {
        setInteractionTable(data){
            this.$refs.interactionTable.setData(data)
        },
        save(){

        }
    },
    mounted(){
        //eventbus控制loading
        let loadingArray = [];
        eventBus.$on('CellChatViewRefreshingStart',()=>{
            loadingArray.push(Loading.service({
                target:".cc-view-container",
                lock:true,
                text:"Refreshing",
                spinner: 'el-icon-loading',
                background: 'rgba(255, 255, 255, 0.8)',
            }));
        })
        eventBus.$on('CellChatViewRefreshingClose',()=>{
            for(let loading of loadingArray)
                loading.close();
        })
    }
};
</script>

<style scoped lang="less">
.cell-chat-container {
    // background-color: white;
    display: flex;
    flex-direction: column;;
    position: relative;
    width: 100%;
    height: 100%;
    // border: 2px solid rgb(200, 200, 200);
    // border-radius: 10px;
    .cc-view-header{
        padding-left:5px;
        flex:0 0 40px;
        display: flex;
        align-items: center;
        background-color: #24292f;
        border:2px solid #24292f;
        #cell-chat-title{
            color:white;
            font-size:20px;
        }   
    }
    .cc-view-container{
        position:relative;
        width: calc(100% - 10px);
        height: 88%;
        display: flex;
        align-items: stretch;
        .cell-chat-view-container{
            flex:1 1 0;
            .cell-chat-view{
                height: 100%;
                width:100%
            }
        }
        .interaction-table-container{
            flex:0 0 240px;
            display: flex;
            align-items: stretch;
            
        }
        .CellChatMode-toggle{
                position:absolute;
                right:0px;
                top:10px;
                /deep/ .el-radio-button__inner {
                    height: 27px;
                    padding:7px 12px;
                }
                z-index:99
        }
        .dataMode-toggle{
                position:absolute;
                right:130px;
                top:10px;
                /deep/ .el-radio-button__inner {
                    height: 27px;
                    padding:7px 12px;
                }
                z-index:99
        }
    }

}


</style>
