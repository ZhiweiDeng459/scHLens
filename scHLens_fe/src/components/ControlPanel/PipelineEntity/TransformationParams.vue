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
                        <el-tooltip
                            content="Normalize each cell by total counts over all genes to let every cell has the same total count(10^6)"
                            placement="right">
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


        <!--RegressOut-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <el-form :label-position="'left'">
                <el-form-item class="form-item">
                    <div style="display:flex;justify-content:space-between;align-items:center;">
                        <div>
                            <b>Regress</b>
                            &nbsp;
                            <el-switch v-model="regressOut" :disabled="mode == 'local'"></el-switch>
                        </div>
                        <el-tooltip
                            content="Take simple linear regression for regressing out 'Total Counts' and 'Mitochondrial cell proportion'"
                            placement="right">
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
                            <b>Scale</b>
                            &nbsp;
                            <el-switch v-model="scale" :disabled="mode == 'local'"></el-switch>
                        </div>
                        <el-tooltip content="scale the expression matrix to unit variance and zero mean" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>

        <div style="padding:5px 5px">
            <i class="el-icon-warning" style="color:orange"></i>
            <a style="word-break: normal; overflow-wrap: normal;">
                scTransform performs both normalization and gene selection. To use scTransform, please select "scTransform"
                as the method in the Gene Selection module. When scTransform is used, scHLens will ignore all other options
                in the Normalization module.
            </a>
        </div>

    </div>
</template>

<script>
import Vue from "vue";
import { Form, FormItem, Input, Tooltip, Checkbox, Card, Switch } from "element-ui";
Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Checkbox.name, Checkbox);
Vue.component(Card.name, Card);
Vue.component(Switch.name, Switch);

export default {
    name: "TransformationParam",
    props: ['mode'],
    data() {
        return {
            log1p: true,
            normalize: true,
            regressOut: false,
            scale: false,
        };
    },
    methods: {
        getParams() {
            let Params = {}

            //log
            if (this.log1p) {
                Params['log1p'] = true;
            }
            //norm
            if (this.normalize) {
                Params['normalize'] = true;
            }
            //regress
            if (this.regressOut) {
                Params['regressOut'] = true;
            }
            //scale
            if (this.scale) {
                Params['scale'] = true;
            }



            return Params;
        }

    },

    watch: {
        /**
         * 监听HighlyVariable的参数变化，使得两类参数不能同时取
         */
        'mode': {
            handler(newValue, oldValue) {
                if (newValue == 'local') {
                    this.regressOut = false;
                    this.scale = false;
                }
            }
        }
    }

};
</script>

<style scoped lang="less">
.form-item {
    margin: 0px;
}</style>
