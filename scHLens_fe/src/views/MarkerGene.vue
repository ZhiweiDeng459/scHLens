<template>
    <div class="marker-gene-container">
        <div class="view-header">
            <div>
                <a id="marker-gene-title">Marker Gene View</a>
                <!-- <el-button type="primary" icon="el-icon-download" style="padding:2px;margin:0px 5px" @click="save"></el-button>
                <el-tooltip style="margin:5px 0px" content="Display the expression of Marker Gene" placement="top">
                    <i class="el-icon-question"></i>
                </el-tooltip> -->
            </div>
        </div>

        <dot-plot ref="dotPlotView" class="dotplot-view" v-show="activeFlag"></dot-plot>
    </div>
</template>

<script>
import DotPlot from "@/components/MarkerGene/DotPlot";
import Vue from "vue";
import { Input, Select, Option, RadioGroup, RadioButton, Loading } from "element-ui";
import eventBus from "@/utils/eventBus.js"


Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(RadioGroup.name, RadioGroup);
Vue.component(RadioButton.name, RadioButton);

export default {
    name: "MarkerGene",
    components: {
        DotPlot,
    },
    data() {
        return {
        };
    },
    computed: {
        curData(){
            return this.$store.state.curData;
        },
        activeFlag(){
            return this.curData.activeFlag['MarkerGene'];
        },
        marker(){
            return this.curData.MK
        }
    },
    methods: {
        save(){
            this.$refs.dotPlotView.saveToFile();
        }
    },
    mounted(){
        //eventbus控制loading
        let loadingArray = [];
        eventBus.$on('MarkerGeneViewRefreshingStart',()=>{
            loadingArray.push(Loading.service({
                target:".dotplot-view",
                lock:true,
                text:"Refreshing",
                spinner: 'el-icon-loading',
                background: 'rgba(255, 255, 255, 0.8)',
            }));
        })
        eventBus.$on('MarkerGeneViewRefreshingClose',()=>{
            for(let loading of loadingArray)
                loading.close();
        })
    }
};
</script>

<style scoped lang="less">
.marker-gene-container{
    width:100%;
    height:100%;
    display: flex;
    // background-color: white;
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
        border:2px solid black;
        #marker-gene-title{
            color:white;
            font-size:20px;
        }
        .marker-gene-panel {
            span {
                user-select: none;
            }
        }          
    }
    .dotplot-view {
        min-height: 0;
        flex: 1 1;
    }
}


/deep/ .el-icon-loading{
    font-size:30px;
}

/deep/ .el-loading-text{
    font-size:25px;
}


</style>
