<template>
    <div>

        
        <!--normalize-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <el-form :label-position="'left'">
                <el-form-item class="form-item">
                    <div style="display:flex;justify-content:space-between;align-items:center;">
                        <div>
                            <b>CPM Normalize</b>
                            &nbsp;
                            <el-switch v-model="normalize"></el-switch>
                        </div>
                        <el-tooltip content="Normalize each cell by total counts over all genes to let every cell has the same total count(10^6)" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>

        <!--log-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <el-form :label-position="'left'">
                <el-form-item class="form-item">
                    <div style="display:flex;justify-content:space-between;align-items:center;">
                        <div>
                            <b>Logarithmize</b>
                            &nbsp;
                            <el-switch v-model="log1p"></el-switch>
                        </div>
                        <el-tooltip content="Take the logarithm for the data" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>

        <!--Gene Selection-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            
            <div slot="header" style="display:flex;flex-direction:column;justify-content:space-between;align-items:center;">
                <div style="border-bottom: 3px solid black;margin-bottom:10px">
                    <b style="font-size:17px">Gene Selection</b>
                </div>
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

        <!--RegressOut-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <el-form :label-position="'left'">
                <el-form-item class="form-item">
                    <div style="display:flex;justify-content:space-between;align-items:center;">
                        <div>
                            <b>Regress</b>
                            &nbsp;
                            <el-switch v-model="regressOut"></el-switch>
                        </div>
                        <el-tooltip content="Take simple linear regression for the data" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>

        <!--Scale-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <el-form :label-position="'left'">
                <el-form-item class="form-item">
                    <div style="display:flex;justify-content:space-between;align-items:center;">
                        <div>
                            <b>scale</b>
                            &nbsp;
                            <el-switch v-model="scale"></el-switch>
                        </div>
                        <el-tooltip content="scale the expression matrix to unit variance and zero mean" placement="right">
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
import { Form, FormItem, Input, Tooltip, Checkbox, Card, Switch} from "element-ui";

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Checkbox.name, Checkbox);
Vue.component(Card.name, Card);
Vue.component(Switch.name, Switch);

export default {
    name: "NormalizeParam",
    props:['mode'],
    data() {
        return {
            log1p:true,
            normalize: true,
            regressOut: false,
            scale: false,
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
            let Params = {}
            //feature selection
            let FSParams = {}
            if(this.feature == 'HighlyVariable'){
                FSParams['highlyVariableGenes'] = {};
                if(this.featureParams['HighlyVariable']['topGenes'] != '')
                    FSParams['highlyVariableGenes']['topGenes'] = Number(this.featureParams['HighlyVariable']['topGenes']);
                if(this.featureParams['HighlyVariable']['minMean'] != '')
                    FSParams['highlyVariableGenes']['minMean'] = Number(this.featureParams['HighlyVariable']['minMean']);
                if(this.featureParams['HighlyVariable']['maxMean'] != '')
                    FSParams['highlyVariableGenes']['maxMean'] = Number(this.featureParams['HighlyVariable']['maxMean']);
                if(this.featureParams['HighlyVariable']['minDisp'] != '')
                    FSParams['highlyVariableGenes']['minDisp'] = Number(this.featureParams['HighlyVariable']['minDisp']);
                if(this.featureParams['HighlyVariable']['maxDisp'] != '')
                    FSParams['highlyVariableGenes']['maxDisp'] = Number(this.featureParams['HighlyVariable']['maxDisp']);
            }
            else if(this.feature == 'scry'){
                FSParams['scry'] = {};
                if(this.featureParams['scry']['topGenes'] != '')
                    FSParams['scry']['topGenes'] = Number(this.featureParams['scry']['topGenes']);
            }
            else if(this.feature == 'SCTransform'){
                FSParams['SCTransform'] = {};
                if(this.featureParams['SCTransform']['topGenes'] != '')
                    FSParams['SCTransform']['topGenes'] = Number(this.featureParams['SCTransform']['topGenes']);
                else{
                    FSParams['SCTransform']['topGenes'] = 1500
                }

            }
            else if(this.feature == 'marker'){
                FSParams['marker'] = {}
            }
            else if(this.feature == 'AllGenes'){
                FSParams['AllGenes'] = {}

            }
            Params['FS'] = FSParams;
            //log
            if(this.log1p){
                Params['log1p'] = true;
            }
            //norm
            if(this.normalize){
                Params['normalize'] = true;
            }
            //regress
            if(this.regressOut){
                Params['regressOut'] = true;
            }
            //scale
            if(this.scale){
                Params['scale'] = true;
            }

            

            return Params;
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
