<template>
    <div>
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="SamplingMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in SamplingMethodOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
            </div >

            <!--Random-->
            <div v-if="SamplingMethod=='Random'">
                <el-form label-width="100px" :label-position="'left'" :model="SamplingMethodParams['Random']">
                    <el-form-item class="form-item" label="Ratio">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="SamplingMethodParams['Random']['sampling_radio']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="Sampling ratio of data set" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="Number">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="SamplingMethodParams['Random']['sampling_num']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="Sampling number of data set" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--Stratified-->
            <div v-if="SamplingMethod=='Stratified'">
                <el-form label-width="100px" :label-position="'left'" :model="SamplingMethodParams['Stratified']">
                    <el-form-item class="form-item" label="Ratio">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="SamplingMethodParams['Stratified']['sampling_radio']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="Sampling ratio of data set" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="Number">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="SamplingMethodParams['Stratified']['sampling_num']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="Sampling number of data set" placement="right">
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
    name: "Sampling",
    data() {
        return {
            SamplingMethod: "Random",
            SamplingMethodOptions: [
                {
                    value: "Random",
                    label: "Random",
                },
                {
                    value: "Stratified",
                    label: "Stratified",
                },
            ],
            SamplingMethodParams:{
                'Random':{
                    sampling_num: '',
                    sampling_radio: '0.2',
                },
                'Stratified':{
                    sampling_num: '',
                    sampling_radio: '0.2',
                },
            },

        };
    },
    methods:{
        getParams(){
            /**
             * 注意要把数字字符串转为数字
             */
            let Params = {};
            if(this.SamplingMethod == 'Random'){
                Params['Random'] = {};
                if(this.SamplingMethodParams['Random']['sampling_num'] != '')
                    Params['Random']['sampling_num'] = Number(this.SamplingMethodParams['Random']['sampling_num']);
                if(this.SamplingMethodParams['Random']['sampling_radio'] != '')
                    Params['Random']['sampling_radio'] = Number(this.SamplingMethodParams['Random']['sampling_radio']);
            }
            if(this.SamplingMethod == 'Stratified'){
                Params['Stratified'] = {};
                if(this.SamplingMethodParams['Stratified']['sampling_num'] != '')
                    Params['Stratified']['sampling_num'] = Number(this.SamplingMethodParams['Stratified']['sampling_num']);
                if(this.SamplingMethodParams['Stratified']['sampling_radio'] != '')
                    Params['Stratified']['sampling_radio'] = Number(this.SamplingMethodParams['Stratified']['sampling_radio']);
            }
            return Params;
        }
    },
    watch:{
        /**
         * 监听参数的变化，因为sampling_radio和sampling_number不能同时取值
         */
        'SamplingMethodParams.Random.sampling_num':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.SamplingMethodParams['Random']['sampling_radio'] = ''
                }
            }
        },
        'SamplingMethodParams.Random.sampling_radio':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.SamplingMethodParams['Random']['sampling_num'] = ''
                }
            }
        },
        'SamplingMethodParams.Stratified.sampling_num':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.SamplingMethodParams['Stratified']['sampling_radio'] = ''
                }
            }
        },
        'SamplingMethodParams.Stratified.sampling_radio':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.SamplingMethodParams['Stratified']['sampling_num'] = ''
                }
            }
        },


    }

};
</script>

<style scoped lang="less">
.form-item {
    margin:0px;
}

</style>
