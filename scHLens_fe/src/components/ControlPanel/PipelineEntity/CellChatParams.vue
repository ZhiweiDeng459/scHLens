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
                    <el-form-item class="form-item" label="Library">
                        <div style="display:flex;justify-content:space-between;align-items:center;">
                            <!-- <el-select v-model="CellChatParams['CellChat']['DatabaseType']" placeholder="Select..." size="mini" style="width:120px">
                                <el-option v-for="item in DatabaseTypeOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                            </el-select> -->
                            <el-cascader
                                size="mini"
                                v-model="CellChatParams.CellChat.curDB"
                                :options="CellChatDBInfo"
                                style="width:180px;margin-left: -40px;margin-right:5px;margin-top:7px">
                                
                            </el-cascader>
                            <el-tooltip content="The type of Receptor-ligand library" placement="right" style="margin-top:7px">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--CellPhoneDB-->
            <div v-if="CellChatMethod=='CellPhoneDB'">
                <el-form label-width="100px" :label-position="'left'">
                    <el-form-item class="form-item" label="Library">
                        <div style="display:flex;justify-content:space-between;align-items:center;">
                            <!-- <el-select v-model="CellChatParams['CellChat']['DatabaseType']" placeholder="Select..." size="mini" style="width:120px">
                                <el-option v-for="item in DatabaseTypeOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                            </el-select> -->
                            <el-cascader
                                size="mini"
                                v-model="CellChatParams.CellPhoneDB.curDB"
                                :options="CellChatDBInfo"
                                style="width:180px;margin-left: -40px;margin-right:5px;margin-top:7px">
                                
                            </el-cascader>
                            <el-tooltip content="The type of Receptor-ligand library" placement="right" style="margin-top:7px">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>

            </div>    

        </el-card>

    </div>
</template>

<script>
import Vue from "vue";
import { Form, FormItem, Input, Select, Option, Radio, Tooltip, Card, Cascader} from "element-ui";
import {queryCellChatDB} from '@/utils/interface.js'

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Card.name, Card);
Vue.component(Cascader.name, Cascader);
export default {
    name: "CellChatParams",
    data() {
        return {
            CellChatDBInfo:[
            ],
            CellChatMethod:'CellChat', //CellChat or CellPhoneDB
            CellChatMethodOptions:[
                {
                    value: "CellChat",
                    label: "CellChat",
                },
                {
                    value: "CellPhoneDB",
                    label: "CellPhoneDB",
                },
            ],
            CellChatParams:{
                'CellChat':{
                    // 'DatabaseType':'human',
                    curDB:[],
                },
                'CellPhoneDB':{
                    curDB:[],
                },
            },
        };
    },
    methods:{
        getParams(){
            let Params = {}
            if(this.CellChatMethod == 'CellChat'){
                Params['CellChat'] = this.CellChatParams['CellChat']
            }
            else if(this.CellChatMethod == 'CellPhoneDB'){
                Params['CellPhoneDB'] = this.CellChatParams['CellPhoneDB']
            }
            return Params
        },

        async updateCellChatDB(){
            await queryCellChatDB().then((res)=>{
                let db_data = res.data

                let CellChatDBInfo = []
                for(let org in db_data){
                    let org_item = {
                        'value':org,
                        'label':org,
                        'children':[]
                    }
                    for(let db_name of db_data[org]){
                        org_item['children'].push({
                            'value':db_name,
                            'label':db_name
                        })
                    }

                    CellChatDBInfo.push(org_item)

                }
                this.CellChatDBInfo = CellChatDBInfo


                //把CellChatDBInfo的第一个值设为初始值 TODO:获取的数据库部分org为空的情况
                for(let method in this.CellChatParams){
                    if(this.CellChatDBInfo.length>0){
                        this.CellChatParams[method]['curDB'] = [this.CellChatDBInfo[0]['value'],this.CellChatDBInfo[0]['children'][0]['value']]
                    }
                }

            }).catch((err)=>{
                console.log('err',err)
            })
        }
    },
    mounted(){
        //更新CellChatDB
        this.updateCellChatDB();
    }
};
</script>

<style scoped lang="less">

.form-item {
    margin:0px;
}

</style>
