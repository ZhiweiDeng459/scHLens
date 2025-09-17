<template>
    <div>

        <!--Gene Selection-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            
            <div slot="header" style="display:flex;flex-direction:column;justify-content:space-between;align-items:center;">
                <div style="display:flex;justify-content:space-between;align-items:center;width:100%">
                    <b>Method</b>
                    <el-select v-model="feature" placeholder="Select..." size="mini" style="width:150px">
                        <el-option v-for="item in featureOptions" :key="item.value" :label="item.label" :value="item.value" :disabled="item.disabled"></el-option>
                    </el-select>
                </div>
            </div>

            <!--AllGenes-->
            <div v-if="feature=='AllGenes'">
                <el-form label-width="100px" :label-position="'left'" :model="featureParams['AllGenes']">
                </el-form>
            </div>   

            <!--HighlyVariable-->
            <div v-if="feature=='HighlyVariable'">
                <el-form label-width="100px" :label-position="'left'" :model="featureParams['HighlyVariable']">
                    <el-form-item class="form-item" label="top Genes">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['HighlyVariable'].topGenes" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The number of top genes (ranked by variability) to keep. This parameter cannot be set together with 'min Mean', 'max Mean', 'min Dispersion', or 'max Dispersion'." placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="min Mean">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['HighlyVariable'].minMean" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The lower bound on Mean for a gene. This parameter cannot be set together with 'top Genes'" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="max Mean">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['HighlyVariable'].maxMean" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The upper bound on Mean for a gene. This parameter cannot be set together with 'top Genes'" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="min Disp">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['HighlyVariable'].minDisp" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The lower bound on Dispersion for a gene. This parameter cannot be set together with 'top Genes'" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="max Disp">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['HighlyVariable'].maxDisp" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The upper bound on Dispersion for a gene. This parameter cannot be set together with 'top Genes'" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
                <div style="padding:5px 5px">
                        <i class="el-icon-warning" style="color:orange"></i>
                        <a style="word-break: normal; overflow-wrap: normal;">
                            Executing HighlyVariable requires prior "Logarithmize"; otherwise, an error may occur.
                        </a>
                </div>
            </div>

            <!--Scry-->
            <div v-if="feature=='scry'">
                <el-form label-width="100px" :label-position="'left'" :model="featureParams['scry']">
                    <el-form-item class="form-item" label="top Genes">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['scry'].topGenes" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="Number of highly-value genes to keep" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--SCTransform-->
            <div v-if="feature=='SCTransform'">
                <el-form label-width="100px" :label-position="'left'" :model="featureParams['SCTransform']">
                    <el-form-item class="form-item" label="top Genes">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="featureParams['SCTransform'].topGenes" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The number of genes as variable features after ranking by residual variance." placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <div style="padding:5px 5px">
                        <i class="el-icon-warning" style="color:orange"></i>
                        <a style="word-break: normal; overflow-wrap: normal;">
                            Since scTransform includes the processes of normalization and gene selection, scHLens will ignore other options in the Normalization module when using scTransform.
                        </a>
                    </div>
                    <div style="padding:5px 5px">
                        <i class="el-icon-warning" style="color:orange"></i>
                        <a style="word-break: normal; overflow-wrap: normal;">
                            scTransform may take a very long time while dealing with large datasets.
                        </a>
                    </div>
                </el-form>
            </div>   

            <!--marker-->
            <div v-if="feature=='marker'">
                <el-form label-width="100px" :label-position="'left'" :model="featureParams['SCTransform']">
                    <i></i>
                    <p style="padding:5px;font-family:YaHei;font-size:12px;color:#434343">
                        This method only select the marker genes as the
                    </p>
                </el-form>
            </div>   


        </el-card>

    </div>
</template>

<script>
import Vue from "vue";
import { Form, FormItem, Input, Tooltip, Checkbox, Card, Switch} from "element-ui";
import {getPipelineParamsErrorT} from "@/utils/objectTemplate";
Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Checkbox.name, Checkbox);
Vue.component(Card.name, Card);
Vue.component(Switch.name, Switch);

export default {
    name: "GeneSelectionParams",
    data() {
        return {
            feature: "HighlyVariable",
            featureOptions: [
                {
                    value: "AllGenes",
                    label: "All Genes",
                },
                {
                    value: "HighlyVariable",
                    label: "HighlyVariable",
                },
                {
                    value: "scry",
                    label: "scry",
                },
                {
                    value: "SCTransform",
                    label: "SCTransform",
                },
                // {
                //     value: "marker",
                //     label: "marker",
                //     disabled:this.mode != 'local',
                // },
            ],
            featureParams:{
                'HighlyVariable':{
                    topGenes: '1500',
                    minMean: '',
                    maxMean: '',
                    minDisp: '',
                    maxDisp: '',
                },
                'scry':{
                    topGenes: '1500',
                },
                'SCTransform':{
                    topGenes: '1500',
                },
                'marker':{

                },
                'AllGenes':{

                },
            },
        };
    },
    methods:{
        getParams(){

            let errMessage = getPipelineParamsErrorT();
            let FSParams = {}
            if(this.feature == 'HighlyVariable'){
                FSParams['highlyVariableGenes'] = {};
                //top Genes
                if(this.featureParams['HighlyVariable']['topGenes'] != ''){
                    let num = Number(this.featureParams['HighlyVariable']['topGenes']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - HighlyVariable - top Genes';
                        errMessage['message'] = '"top Genes" should be a valid number';
                        return errMessage;
                    }
                    //判断是否为正整数
                    if(num <= 0 || !Number.isInteger(num)){
                        errMessage['location'] = 'Gene Selection - HighlyVariable - top Genes';
                        errMessage['message'] = '"top Genes" should be a positive integer';
                        return errMessage;
                    }
                    FSParams['highlyVariableGenes']['topGenes'] = num
                }
                //min Mean
                if(this.featureParams['HighlyVariable']['minMean'] != ''){
                    let num = Number(this.featureParams['HighlyVariable']['minMean']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - HighlyVariable - min Mean';
                        errMessage['message'] = '"min Mean" should be a valid number';
                        return errMessage;
                    }
                    FSParams['highlyVariableGenes']['minMean'] = num
                }
                //max Mean
                if(this.featureParams['HighlyVariable']['maxMean'] != ''){
                    let num = Number(this.featureParams['HighlyVariable']['maxMean']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - HighlyVariable - max Mean';
                        errMessage['message'] = '"max Mean" should be a valid number';
                        return errMessage;
                    }
                    FSParams['highlyVariableGenes']['maxMean'] = num
                }
                //min Disp
                if(this.featureParams['HighlyVariable']['minDisp'] != ''){
                    let num = Number(this.featureParams['HighlyVariable']['minDisp']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - HighlyVariable - min Disp';
                        errMessage['message'] = '"min Disp" should be a valid number';
                        return errMessage;
                    }
                    FSParams['highlyVariableGenes']['minDisp'] = num
                }
                //max Disp
                if(this.featureParams['HighlyVariable']['maxDisp'] != ''){
                    let num = Number(this.featureParams['HighlyVariable']['maxDisp']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - HighlyVariable - max Disp';
                        errMessage['message'] = '"max Disp" should be a valid number';
                        return errMessage;
                    }
                    FSParams['highlyVariableGenes']['maxDisp'] = num
                }
                //判断是否一个参数都没有设置
                if(FSParams['highlyVariableGenes']['topGenes'] == undefined && FSParams['highlyVariableGenes']['minMean'] == undefined && FSParams['highlyVariableGenes']['maxMean'] == undefined && FSParams['highlyVariableGenes']['minDisp'] == undefined && FSParams['highlyVariableGenes']['maxDisp'] == undefined){
                    errMessage['location'] = 'Gene Selection - HighlyVariable';
                    errMessage['message'] = 'Please set at least one parameter for "HighlyVariable"';
                    return errMessage;
                }
                
            }

            else if(this.feature == 'scry'){
                FSParams['scry'] = {};
                //top Genes
                if(this.featureParams['scry']['topGenes'] != ''){
                    let num = Number(this.featureParams['scry']['topGenes']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - scry - top Genes';
                        errMessage['message'] = '"top Genes" should be a valid number';
                        return errMessage;
                    }
                    //判断是否为正整数
                    if(num <= 0 || !Number.isInteger(num)){
                        errMessage['location'] = 'Gene Selection - scry - top Genes';
                        errMessage['message'] = '"top Genes" should be a positive integer';
                        return errMessage;
                    }
                    FSParams['scry']['topGenes'] = num
                }
                else{//至少需要设置一个参数
                    errMessage['location'] = 'Gene Selection - scry - top Genes';
                    errMessage['message'] = '"top Genes" should be set';
                    return errMessage;
                }
            }

            else if(this.feature == 'SCTransform'){
                FSParams['SCTransform'] = {};

                //top Genes
                if(this.featureParams['SCTransform']['topGenes'] != ''){
                    let num = Number(this.featureParams['SCTransform']['topGenes']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Gene Selection - SCTransform - top Genes';
                        errMessage['message'] = '"top Genes" should be a valid number';
                        return errMessage;
                    }
                    //判断是否为正整数
                    if(num <= 0 || !Number.isInteger(num)){
                        errMessage['location'] = 'Gene Selection - SCTransform - top Genes';
                        errMessage['message'] = '"top Genes" should be a positive integer';
                        return errMessage;
                    }
                    FSParams['SCTransform']['topGenes'] = num
                }
                else{//至少需要设置一个参数
                    errMessage['location'] = 'Gene Selection - SCTransform - top Genes';
                    errMessage['message'] = '"top Genes" should be set';
                    return errMessage;
                }

            }
            else if(this.feature == 'marker'){
                FSParams['marker'] = {}
            }
            else if(this.feature == 'AllGenes'){
                FSParams['AllGenes'] = {}

            }
            

            return FSParams;
        }

    },

    watch:{
        /**
         * 监听HighlyVariable的参数变化，使得两类参数不能同时取
         */
        'featureParams.HighlyVariable.topGenes':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.featureParams['HighlyVariable']['minMean'] = ''
                    this.featureParams['HighlyVariable']['maxMean'] = ''
                    this.featureParams['HighlyVariable']['minDisp'] = ''
                    this.featureParams['HighlyVariable']['maxDisp'] = ''
                }
            }
        },
        'featureParams.HighlyVariable.minMean':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.featureParams['HighlyVariable']['topGenes'] = ''
                }
            }
        },
        'featureParams.HighlyVariable.maxMean':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.featureParams['HighlyVariable']['topGenes'] = ''
                }
            }
        },
        'featureParams.HighlyVariable.minDisp':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.featureParams['HighlyVariable']['topGenes'] = ''
                }
            }
        },
        'featureParams.HighlyVariable.maxDisp':{
            handler(newValue,oldValue){
                if(newValue!='' && newValue !== undefined && newValue !== null){
                    this.featureParams['HighlyVariable']['topGenes'] = ''
                }
            }
        },
    }

};
</script>

<style scoped lang="less">
.form-item{
    margin:0px;
}
</style>
