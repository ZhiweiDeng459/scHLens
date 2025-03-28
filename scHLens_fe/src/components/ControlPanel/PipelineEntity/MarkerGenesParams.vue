<template>
    <div>
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <el-form label-width="100px" :label-position="'left'">
                <el-form-item class="form-item" label="Method">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-select v-model="markerMethod" placeholder="Select..." size="mini" style="width:110px">
                            <el-option v-for="item in markerOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                        </el-select>
                        <el-tooltip content="DEG Identification method" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
                <el-form-item class="form-item" label="nGenes">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="nGenes" size="mini" style="width: 110px;"></el-input>
                        <el-tooltip content="The number of DEGs" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
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
    name: "MarkerGenesParams",
    data() {
        return {
            markerMethod: "t-test",
            nGenes: '50',
            markerOptions: [
                {
                    value: "t-test",
                    label: "t-test",
                },
                {
                    value: "t-test_overestim_var",
                    label: "t-test(overestimate variance)",
                },
                {
                    value: "wilcoxon-test",
                    label: "wilcoxon-test",
                },
                {
                    value: "wilcoxon-test(TIE)",
                    label: "wilcoxon-test(TIE)",
                },
                {
                    value: "logreg",
                    label: "logistic regression",
                },
                {
                    value: "COSG",
                    label: "COSG",
                },
            ],
        };
    },
    methods:{
        getParams(){
            /**
             * 注意要把数字字符串转为数字
             */
            let Params = {}
            if(this.markerMethod != ''){
                Params['markerMethod'] = this.markerMethod;

            }
            if(this.nGenes != ''){
                Params['nGenes'] = Number(this.nGenes);
            }
            return Params;
        }
    }
};
</script>

<style scoped lang="less">
.form-item {
    margin: 0px;
}

</style>
