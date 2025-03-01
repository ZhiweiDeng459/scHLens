<template>
    <div class="trajectory-inference-container">
        <div class="view-header">
            <div>
                <b id="trajectory-inference-title">Trajectory-Inference View</b>
                <el-button type="primary" icon="el-icon-download" style="padding:2px;margin:0px 4px" @click="save"></el-button>
                <el-tooltip style="margin:5px 2px" content="Display the distruibution of cell-clusters" placement="top">
                    <i class="el-icon-question"></i>
                </el-tooltip>
            </div>
        </div>
        <div id="TrajectoryInferenceViewContainer" v-show="activeFlag">
            <el-radio-group v-show="curMethod=='paga'" v-model="pagaMode" size="mini" class="mode-toggle">
                <el-radio-button label="Plot"></el-radio-button>
                <el-radio-button label="Scatter"></el-radio-button>
            </el-radio-group>
            <el-radio-group v-show="curMethod=='slingshot'" v-model="slingshotMode" size="mini" class="mode-toggle">
                <el-radio-button label="Curve"></el-radio-button>
                <el-radio-button label="Link"></el-radio-button>
            </el-radio-group>
            <el-radio-group v-show="curMethod=='slingshot'" v-model="colorMode" size="mini" class="mode-toggle-2">
                <el-radio-button label="group"></el-radio-button>
                <el-radio-button label="pseudotime"></el-radio-button>
            </el-radio-group>
            <TrajectoryInferencePlot v-if="pagaMode=='Plot' && curMethod == 'paga'" ref="TrajectoryInferencePlot" class="TrajectoryInferenceView"></TrajectoryInferencePlot>
            <TrajectoryInferenceScatter  v-if="pagaMode=='Scatter' && curMethod == 'paga'" ref="TrajectoryInferenceScatter" class="TrajectoryInferenceView"></TrajectoryInferenceScatter>
            <TrajectoryInferenceCurve v-if="slingshotMode=='Curve' && curMethod == 'slingshot'" ref="TrajectoryInferenceCurve" class="TrajectoryInferenceView" :colorMode="colorMode"></TrajectoryInferenceCurve>
            <TrajectoryInferenceLink v-if="slingshotMode=='Link' && curMethod == 'slingshot'" ref="TrajectoryInferenceLink" class="TrajectoryInferenceView" :colorMode="colorMode"></TrajectoryInferenceLink>
        </div>
    </div>
</template>

<script>
import TrajectoryInferencePlot from "@/components/TrajectoryInference/TrajectoryInferencePlot"
import TrajectoryInferenceScatter from "@/components/TrajectoryInference/TrajectoryInferenceScatter"
import TrajectoryInferenceCurve from "@/components/TrajectoryInference/TrajectoryInferenceCurve"
import TrajectoryInferenceLink from "@/components/TrajectoryInference/TrajectoryInferenceLink"


export default {
    name:'TrajectoryInference',
    components:{
        TrajectoryInferencePlot,
        TrajectoryInferenceScatter,
        TrajectoryInferenceCurve,
        TrajectoryInferenceLink,
    },
    data() {
        return {
            pagaMode:'Plot', //Plot or Scatter
            slingshotMode:'Curve', //Curve or Link
            colorMode: 'group'//group or pseudotime // Only for slingshot 
        }
    },
    computed:{
        activeFlag(){
            return this.$store.state.curData.activeFlag['TrajectoryInference'];
        },
        curMethod(){
            if(this.$store.state.dataList.length == 0){//还没有运行过数据
                return undefined;
            }
            else if(!('TI' in this.$store.state.curData.paramsObj)){//当前节点没有执行过轨迹推断
                return undefined;
            }
            else{
                return Object.keys(this.$store.state.curData.paramsObj.TI)[0]
            }
        },

    },
    watch:{

    },
    methods:{
        save(){
            if(this.mode == 'Plot')
                this.$refs.TrajectoryInferencePlot.saveToFile();
            else if(this.mode == 'Scatter')
                this.$refs.TrajectoryInferenceScatter.saveToFile();
        }
    }
}
</script>

<style lang="less">

.trajectory-inference-container{
    width:100%;
    height:100%;
    display: flex;
    // background-color: white;
    flex-direction: column;
    align-items: stretch;
    // border: 2px solid rgb(200, 200, 200);
    // border-radius: 10px;
    .view-header{
        padding:5px;
        border-bottom: 2px solid lightgray;
        display: flex;
        #trajectory-inference-title{
            font-family:YaHei;
        }          
    }
    #TrajectoryInferenceViewContainer{
        position:relative;
        width: calc(100% - 10px);
        height: 92%;
        .TrajectoryInferenceView{
            position:absolute;
            margin: 5px;
            height: 100%;
            width:100%
        }
        .mode-toggle {
            position: absolute;
            z-index: 99;
            right: 0%;
            top: 8px;
            /deep/ .el-radio-button__inner {
                padding: 5px 8px;
            }
        }
        .mode-toggle-2 {
            position: absolute;
            z-index: 99;
            right: 0%;
            top: 40px;
            /deep/ .el-radio-button__inner {
                padding: 5px 8px;
            }
        }
    }

}




</style>