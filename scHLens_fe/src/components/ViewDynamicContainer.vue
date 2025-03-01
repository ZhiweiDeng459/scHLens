<template>
    <div class='outer-container' ref="outer-container">
        <CellProjection v-if="chooseFlags[rank]['CellProjection']"/>
        <GeneProjection v-if="chooseFlags[rank]['GeneProjection']"/>
        <GeneExpression v-if="chooseFlags[rank]['GeneExpression']"/>
        <MarkerGene v-if="chooseFlags[rank]['MarkerGene']"/>
        <TrajectoryInference v-if="chooseFlags[rank]['TrajectoryInference']"/>
        <CellChat v-if="chooseFlags[rank]['CellChat']"/>
        <!--遮罩层-->
        <div v-show="isEmpty()" class="empty-modal-container">
            <el-empty 
                image="icons/view-empty.svg"
                description="NO VIEW" 
                :image-size="200" 
                style="height:100%;width:100%">
            </el-empty>
        </div>
    </div>
</template>

<script>
import Vue from "vue";
import CellProjection from "@/views/CellProjection";
import GeneProjection from "@/views/GeneProjection";
import GeneExpression from "@/views/GeneExpression";
import MarkerGene from "@/views/MarkerGene";
import TrajectoryInference from  "@/views/TrajectoryInference"
import CellChat from "@/views/CellChat"

//动态的视图容器，可以根据chooseFlag切换
export default {
    name: "ViewDynamicContainer",
    props:['rank'],
    components: {
        CellProjection,
        GeneProjection,
        GeneExpression,
        MarkerGene,
        TrajectoryInference,
        CellChat,
    },
    data(){
        return{
        }
    },
    computed:{
        chooseFlags(){
            return this.$store.state.chooseFlags;
        },
        curData(){
            return this.$store.state.curData;
        }

    },
    methods:{
        isEmpty(){//判断容器是否为空
            if(this.rank != 0 && this.rank != 1)
                return true;
            for(let key in this.chooseFlags[this.rank]){
                
                if(this.chooseFlags[this.rank][key]){
                    return false
                }
            }
            return true
        },
    },
    watch:{
        curData(){
            //保持长宽比
            let height = this.$refs['outer-container'].clientHeight;
            this.$refs['outer-container'].style.width=String(height) + "px"
        }
    }
}
</script>

<style style scoped lang="less">
    .outer-container{
        background-color: white;
        border: 2px solid rgb(200, 200, 200);
        border-radius: 2px;
        height: 100%;
        min-height: 0;
        min-width: 0;
        // height:600px;
        // min-height: 600px;
        // width: 600px;
        // min-width: 600px;
        // aspect-ratio: 1;//固定长宽比
    }
    .empty-modal-container{
        height: 100%;
        width: 100%;
        /deep/ .el-empty__description p{
            font-size:22px
        }
    }
</style>