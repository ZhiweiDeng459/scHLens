<template>
    <div class="gene-projection-container">
        <div class="view-header">
            <div>
                <a id="gene-projection-title">Gene Projection View</a>
            </div>
        </div>

        <div class="view-container">
            <el-radio-group v-model="mode" size="mini" class="mode-toggle" v-show="activeFlag">
                <el-radio-button label="Scatter"></el-radio-button>
                <el-radio-button label="Density"></el-radio-button>
            </el-radio-group>
            <gene-scatter ref="geneScatter" class="gene-scatter-view" v-if="activeFlag&&mode=='Scatter'"></gene-scatter>
            <density-scatter ref="densityScatter" class="density-scatter-view" v-if="activeFlag&&mode=='Density'"></density-scatter>
        </div>
    </div>
</template>

<script>
import GeneScatter from "@/components/GeneProjection/GeneScatter";
import DensityScatter from "@/components/GeneProjection/DensityScatter";
import Vue from "vue";
import { Input, Button, RadioGroup, RadioButton,Loading} from "element-ui";
import eventBus from "@/utils/eventBus.js"


Vue.component(Input.name, Input);
Vue.component(Button.name, Button);
Vue.component(RadioGroup.name, RadioGroup);
Vue.component(RadioButton.name, RadioButton);

export default {
    name: "GeneProjection",
    components: {
        GeneScatter,
        DensityScatter,
    },
    computed:{
        curData(){
            return this.$store.state.curData;
        },
        cellData(){
            return this.curData.cellData;
        },
        dataset(){
            return this.curData.paramsObj['dataset'];
        },
        activeFlag(){
            return this.curData.activeFlag['GeneProjection'];
        },
        JobId(){
            return this.$store.state.JobId;
        }

    },
    data() {
        return {
            searchGeneName: "",
            mode: "Scatter",
        };
    },
    methods: {
        save(){
            if(this.mode == 'Scatter')
                this.$refs.geneScatter.saveToFile();
            else if(this.mode == 'Density')
                this.$refs.densityScatter.saveToFile();
        }
    },
    mounted(){
        //eventbus控制loading
        let loadingArray = [];
        eventBus.$on('GeneProjectionViewRefreshingStart',()=>{
            loadingArray.push(Loading.service({
                target:".view-container",
                lock:true,
                text:"Refreshing",
                spinner: 'el-icon-loading',
                background: 'rgba(255, 255, 255, 0.8)',
            }));
        })
        eventBus.$on('GeneProjectionViewRefreshingClose',()=>{
            for(let loading of loadingArray)
                loading.close();
        })
    }
};
</script>

<style scoped lang="less">
.gene-projection-container{
    // background-color: white;
    position: relative;
    width: 100%;
    height: 100%;
    display: flex;
    flex-direction: column;
    // border: 2px solid rgb(200, 200, 200);
    // border-radius: 10px;
    .view-header{
        padding-left:5px;
        flex:0 0 40px;
        display: flex;
        align-items: center;
        background-color: #24292f;
        border:2px solid #24292f;
        #gene-projection-title{
            color:white;
            font-size:20px;
        }   
    }
    .gene-panel {
            display: flex;
            align-items: center;
            /deep/ .el-input__inner {
                height: 22px;
            }
        .recommendGenePopoverInGP{
            height:20px;
            margin-left: 2px;
        }
    }
    .view-container{
        position: relative;
        flex: 1 1 0;
        .gene-scatter-view {
            // position:absolute;
            width: 100%;
            height: 100%;
        }
        .density-scatter-view {
            // position:absolute;
            width: 100%;
            height: 100%;
        }
        .mode-toggle {
            position: absolute;
            z-index: 99;
            right: 12%;
            top: 8px;
            /deep/ .el-radio-button__inner {
                padding: 5px 8px;
            }
        }
        .currentGeneLabel{
            position:absolute;
            z-index:98;
            top:3px;
            left:10px
        }
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
