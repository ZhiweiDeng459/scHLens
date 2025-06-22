<template>
    <div class='outer-container' ref="outer-container">
        <CellProjection v-if="!isEmpty()&&type=='CellProjection'"/>
        <GeneProjection v-if="!isEmpty()&&type=='GeneProjection'"/>
        <GeneExpression v-if="!isEmpty()&&type=='GeneExpression'"/>
        <MarkerGene v-if="!isEmpty()&&type=='MarkerGene'"/>
        <TrajectoryInference v-if="!isEmpty()&&type=='TrajectoryInference'"/>
        <CellChat v-if="!isEmpty()&&type=='CellChat'"/>
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

//静态的视图容器，视图种类只能在一开始设定，不能切换
export default {
    name: "ViewContainer",
    props:['type'],
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
        },
        marker(){
            return this.curData.MK;
        },

    },
    methods:{
        isEmpty(){//判断是否有内容
            if(this.type == 'MarkerGene'){
                if(this.marker === undefined || this.marker === null || this.marker.length == 0){
                    return true
                }
            }
            return false

        },
    },
    watch:{
    }
}
</script>

<style style scoped lang="less">
    .outer-container{
        background-color: white;
        border: 2px solid rgb(200, 200, 200);
        border-radius: 2px;
        // height: 100%;
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