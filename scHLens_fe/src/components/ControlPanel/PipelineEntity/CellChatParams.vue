<template>
    <div>
        <!-- <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Database Type</b>
                <el-select v-model="DatabaseType" placeholder="Select..." size="mini" style="width:130px">
                    <el-option v-for="item in DatabaseTypeOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
            </div>
        </el-card> -->

        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="CellChatMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in CellChatMethodOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
            </div >

            <!--CellChat-->
            <div v-if="CellChatMethod=='CellChat'">
                <el-form label-width="100px" :label-position="'left'">
                    <el-form-item class="form-item" label="database">
                        <div style="display:flex;justify-content:space-between;align-items:center;">
                            <el-select v-model="CellChatParams['CellChat']['DatabaseType']" placeholder="Select..." size="mini" style="width:150px">
                                <el-option v-for="item in DatabaseTypeOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                            </el-select>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--NicheNet-->
            <div v-if="CellChatMethod=='NicheNet'">
                
            </div>    

        </el-card>

    </div>
</template>

<script>
import Vue from "vue";
import { Form, FormItem, Input, Select, Option, Radio, Tooltip, Card} from "element-ui";

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Card.name, Card);

export default {
    name: "CellChatParams",
    data() {
        return {
            CellChatMethod:'CellChat', //CellChat or NicheNet
            CellChatMethodOptions:[
                {
                    value: "CellChat",
                    label: "CellChat",
                },
                {
                    value: "NicheNet",
                    label: "NicheNet",
                },
            ],
            CellChatParams:{
                'CellChat':{
                    'DatabaseType':'human',
                },
                'NicheNet':{

                },
            },

            DatabaseTypeOptions: [
                {
                    value: "human",
                    label: "human",
                },
                {
                    value: "mouse",
                    label: "mouse",
                },
            ],
        };
    },
    methods:{
        getParams(){
            let Params = {}
            if(this.CellChatMethod == 'CellChat'){
                Params['CellChat'] = this.CellChatParams['CellChat']
            }
            else if(this.CellChatMethod == 'NicheNet'){
                Params['NicheNet'] = this.CellChatParams['NichNet']
            }
            return Params
        }
    }
};
</script>

<style scoped lang="less">

.form-item {
    margin:0px;
}

</style>
