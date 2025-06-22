<template>
    <div>
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="TIMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in TIMethodOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
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
    name: "TrajectoryInferenceParams",
    data() {
        return {
            TIMethod:'paga',
            TIMethodOptions: [
                {
                    value: "paga",
                    label: "paga",
                },
                {
                    value: "slingshot",
                    label: "slingshot",
                },
            ],
        };
    },
    methods:{
        getParams(){
            let Params = {}
            if(this.TIMethod == 'paga'){
                Params['paga'] = {}
            }
            else if(this.TIMethod == 'slingshot'){
                Params['slingshot'] = {}
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
