<template>

    <div v-show="modal_visible" class="modal">

        <div class="info">
            <div class="top-info">
                <b>{{ top_info }}</b>
            </div>

            <div class="progress">
                <el-progress class="progress-bar" :percentage="percentage" :stroke-width="24" :show-text="false"></el-progress>
                <b class="progress-text">{{ `${Math.round(percentage)}%` }}</b>
            </div>


            <div class="bottom-info">
                <b>{{bottom_info}}</b>
            </div>
        </div>

    </div>

</template>

<script>
import Vue from "vue";

import {Progress} from "element-ui";
Vue.component(Progress.name, Progress);


export default {
    name: "LoadingModal",
    data(){
        return {
            modal_visible:false,
            top_info:'',
            bottom_info:'',
            percentage:50,
        }
    },
    methods:{
        /**
         * 外部接口
         */
        setVisible(modal_visible){
            this.modal_visible = modal_visible
        },
        setTopInfo(text){
            this.top_info = text
        },
        setBottomInfo(text){
            this.bottom_info = text
        },
        setPercentage(percentage){
            this.percentage = percentage
        }


        /**
         * 内部函数
         */
        

    }
}
</script>

<style scoped lang="less">
    .modal{
        background-color:white;
        width: 100vw;
        height: 100vh;
        opacity:0.95;
        z-index: 999;
        position: fixed;
        top:0px;
        left: 0px;
        display: flex;
        flex-direction: column;
        align-items: center;

        .info{
            width:600px;
            display: flex;
            margin-top:18%;
            flex-direction: column;
            align-items: stretch;
            justify-content: center;
            .top-info{
                margin-bottom: 15px;
                font-size:25px;
                // font-family: 'YaHei';
            }

            .progress{
                display: flex;
                flex-direction: row;
                align-items: center;

                .progress-bar{
                    width: 530px;
                }

                .progress-text{
                    font-size: 20px;
                    margin-left: 10px;
                }
            }

            .bottom-info{
                margin-top: 8px;
                font-size:20px;
                // font-family: 'YaHei';

            }
        }
    }


</style>