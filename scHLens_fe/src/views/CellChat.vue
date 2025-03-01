<template>
        <div class="cell-chat-container">
            <div class="view-header">
                <div style="flex:1 1 0">
                    <b id="cell-chat-title">Cell Chat View</b>
                    <el-button type="primary" icon="el-icon-download" style="padding:2px;margin:0px 5px" @click="save"></el-button>
                    <el-tooltip content="Display the difference of gene-expression among clusters" placement="top">
                        <i class="el-icon-question"></i>
                    </el-tooltip>
                </div>
            </div>
            <div class="cell-chat-view-container">
                <el-radio-group v-show="curMethod=='CellChat'" v-model="CellChatMode" size="mini" class="mode-toggle">
                    <el-radio-button label="Chord"></el-radio-button>
                    <el-radio-button label="Edge"></el-radio-button>
                </el-radio-group>
                <CellChatChordPlot v-if="curMethod=='CellChat' && CellChatMode=='Chord'" class="cell-chat-view"></CellChatChordPlot>
                <CellChatEdgePlot v-if="curMethod=='CellChat' && CellChatMode=='Edge'"  class="cell-chat-view"></CellChatEdgePlot>
            </div>
        </div>
</template>

<script>
import CellChatChordPlot from "@/components/CellChat/CellChatChordPlot"
import CellChatEdgePlot from "@/components/CellChat/CellChatEdgePlot"
export default {
    name: "CellChat",
    components: {
        CellChatChordPlot,
        CellChatEdgePlot
    },
    data(){
        return {
            CellChatMode:'Chord' //Chord or Edge
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

        save(){

        }
    },
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
    .view-header{
        padding:5px;
        border-bottom: 2px solid lightgray;
        display: flex;
        #cell-chat-title{
            font-family:YaHei;
        }   
    }
    .cell-chat-view-container{
        position:relative;
        width: calc(100% - 10px);
        height: 92%;
        .cell-chat-view{
            position:absolute;
            margin: 5px;
            height: 100%;
            width:100%
        }
        .mode-toggle{
                position:absolute;
                right:10px;
                top:10px;
                /deep/ .el-radio-button__inner {
                    height: 27px;
                }
                z-index:99
        }
    }

}


</style>
