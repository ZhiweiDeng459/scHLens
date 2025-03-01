<template>
    <!-- 图形选择栏 -->
    <div class='choose-container'>
        <el-checkbox-group size="small" :max="2" v-model="chooseResult" >
            <el-checkbox-button
                
                v-for="view in views" 
                @change="allocateContainer(view['value'])" 
                :disabled="!activeFlag[view['value']]"
                :label="view['value']" 
                :key="view['value']"
                >{{view['label']}}</el-checkbox-button>
        </el-checkbox-group>


    </div>  
</template>

<script>
import Vue from "vue";
import { CheckboxGroup,CheckboxButton} from "element-ui";
Vue.use(CheckboxGroup)
Vue.use(CheckboxButton)
export default {
    name: "ViewChooser",
    computed:{
        activeFlag(){
            return this.$store.state.curData.activeFlag;
        },
        chooseFlags(){
            return this.$store.state.chooseFlags;
        },
    },
    methods:{
        isIdle(index){
            //判断index号容器是否空闲
            for(let key in this.chooseFlags[index]){
                if(this.chooseFlags[index][key])
                    return false;
            }
            return true;
        },
        allocateContainer(comp){//分配容器
            //先检测是否使用了容器
            for(let i = 0;i < this.chooseFlags.length;i++){
                if(this.chooseFlags[i][comp]){
                    this.$store.commit("unchooseView",[i, comp]);
                    return ;
                }
            }
            //分配一个容器
            let flag0 = true
            for(let x in this.chooseFlags[0]){
                if(this.chooseFlags[0][x]){
                    flag0 = false;
                    break;
                }
            }
            if(flag0){
                this.$store.commit("chooseView",[0, comp]);
                return ;
            }
            this.$store.commit("chooseView",[1, comp]);

            return ;
        }

    },
    watch:{
        chooseFlags:{
            'deep':true,
            handler(){
                this.chooseResult.length = 0
                for(let i = 0;i < this.chooseFlags.length;i++){
                    for(let key in this.chooseFlags[i]){
                        if(this.chooseFlags[i][key]){
                            this.chooseResult.push(key)
                        }
                    }
                }    
            }
        }
    },
    data(){
        return{
            views:[
                {
                  label:'Cell Projection',
                  value:'CellProjection', 
                },
                {
                  label:'Gene Projection',
                  value:'GeneProjection', 
                },
                {
                  label:'Gene Expression',
                  value:'GeneExpression', 
                },
                // {
                //   label:'Marker Gene',
                //   value:'MarkerGene', 
                // },
                // {
                //   label:'Trajectory Inference',
                //   value:'TrajectoryInference', 
                // },
                // {
                //   label:'Cell Chat',
                //   value:'CellChat', 
                // },
            ],
            chooseResult:[],
        }
    }
}
</script>

<style scoped lang="less">

    .choose-container{
        display: flex;
        align-items: flex-end;
        // margin-left: 5px;
        /deep/ .el-checkbox-button__inner{
            border-color:rgb(200, 200, 200) !important;
            border-bottom: none !important;
            border-radius: 0px !important;
        }
    }
</style>