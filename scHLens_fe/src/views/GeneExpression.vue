<template>
        <div class="gene-expression-container">
            <div class="view-header">
                <div>
                    <a id="gene-expression-title">Distribution of Gene Expression View</a>
                    <!-- <el-button type="primary" icon="el-icon-download" style="padding:2px;margin:0px 5px" @click="save"></el-button>
                    <el-tooltip content="Display the difference of gene-expression among clusters" placement="top">
                        <i class="el-icon-question"></i>
                    </el-tooltip> -->
                </div>
            </div>
            <div             
                class="violin-container"  v-show="activeFlag">
                <el-radio-group v-model="mode" size="mini" class="mode-toggle">
                        <el-radio-button label="Box"></el-radio-button>
                        <el-radio-button label="Strip"></el-radio-button>
                </el-radio-group>
                <violin ref="violinView" class="violin-view" :mode = "mode"></violin>
            </div>
        </div>
</template>

<script>
import Violin from "@/components/GeneExpression/Violin";
import { Form, Tabs, TabPane, Loading, Drawer,Message } from "element-ui";

import eventBus from "@/utils/eventBus.js"

export default {
    name: "GeneExpression",
    props:['data','show'],
    data(){
        return {
            mode:"Box",
        }
    },
    components: {
        Violin,
    },
    computed:{
        curData() {
            return this.$store.state.curData;
        },
        activeFlag(){
            return this.curData.activeFlag['GeneExpression'];
        },
    },
    methods:{
        save(){
            this.$refs.violinView.saveToFile();
        }
    },
    mounted(){
        //eventbus控制loading
        let loadingArray = [];
        eventBus.$on('GeneExpressionViewRefreshingStart',()=>{
            loadingArray.push(Loading.service({
                target:".violin-container",
                lock:true,
                text:"Refreshing",
                spinner: 'el-icon-loading',
                background: 'rgba(255, 255, 255, 0.8)',
            }));
        })
        eventBus.$on('GeneExpressionViewRefreshingClose',()=>{
            for(let loading of loadingArray)
                loading.close();
        })
    }
};
</script>

<style scoped lang="less">
.gene-expression-container{
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
        border:2px solid #24292f;
        #gene-expression-title{
            // font-family:YaHei;
            color:white;
            font-size:20px;
        }  
    }
    .violin-container{
        width: 100%;
        height: 100%;
        position: relative;
        .violin-view{
            width:100%;
            height:100%;
            
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

//调整loading的字体大小

/deep/ .el-icon-loading{
    font-size:30px;
}

/deep/ .el-loading-text{
    font-size:25px;
}

</style>
