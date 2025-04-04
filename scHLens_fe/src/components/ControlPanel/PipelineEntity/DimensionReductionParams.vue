 <template>
    <div>
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="drMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in DRMethodOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
            </div >

            <!--UMAP-->
            <div v-if="drMethod=='UMAP'">
                <el-form label-width="100px" :label-position="'left'" :model="DRParams['UMAP']">
                    <el-form-item class="form-item" label="pre-DR">
                        <div style="display:flex;justify-content:space-between;align-items:center;height:40px">
                            <el-switch v-model="DRParams['UMAP']['PreDR']"></el-switch>
                            <el-tooltip content="Pre-Dimensionality Reduction using PCA" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="min Dist">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="DRParams['UMAP'].minDist" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The minimum distance between points in the two-dimensional embedding" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="nNeighbors">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="DRParams['UMAP'].n_neighbors" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The number of points in the local neighborhood used for manifold learning" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--t-SNE-->
            <div v-if="drMethod=='T-SNE'">
                <el-form label-width="100px" :label-position="'left'" :model="DRParams['T-SNE']">
                    <el-form-item class="form-item" label="pre-DR">
                        <div style="display:flex;justify-content:space-between;align-items:center;height:40px">
                            <el-switch v-model="DRParams['T-SNE']['PreDR']"></el-switch>
                            <el-tooltip content="Pre-Dimensionality Reduction using PCA" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="perplexity">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="DRParams['T-SNE'].perplexity" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The effective number of neighbors considered for each data point when computing pairwise similarities in the high-dimensional space" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--PCA-->
            <div v-if="drMethod=='PCA'">
                <el-form label-width="100px" :label-position="'left'">
                </el-form>
            </div>

        </el-card>
    </div>
</template>

<script>
import Vue from "vue";
import { Form, FormItem, Input, Select, Option, Radio, Tooltip } from "element-ui";
import {getPipelineParamsErrorT} from "@/utils/objectTemplate";

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);

export default {
    name: "DimensionReductionParams",
    data() {
        return {
            drMethod: "UMAP",
            DRParams:{
                'PCA':{

                },
                'UMAP':{
                    minDist:'0.5',
                    n_neighbors:15,
                    PreDR:false

                },
                'T-SNE':{
                    perplexity:'30',
                    PreDR:false
                }
            },
            DRMethodOptions: [
                {
                    value: "PCA",
                    label: "PCA",
                },
                {
                    value: "T-SNE",
                    label: "T-SNE",
                },
                {
                    value: "UMAP",
                    label: "UMAP",
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
            let errMessage = getPipelineParamsErrorT()
            if(this.drMethod == 'PCA'){
                Params['PCA'] = {}
            }
            else if(this.drMethod == 'UMAP'){
                Params['UMAP'] = {}
                //minDist
                if(this.DRParams['UMAP']['minDist'] != ''){
                    let num = Number(this.DRParams['UMAP']['minDist']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Visualization - UMAP - min Dist';
                        errMessage['message'] = '"min Dist" should be a valid number';
                        return errMessage;
                    }
                    //判断是否为正数
                    if(num <= 0){
                        errMessage['location'] = 'Visualization - UMAP - min Dist';
                        errMessage['message'] = '"min Dist" should be a positive number';
                        return errMessage;
                    }
                    Params['UMAP']['minDist'] = num
                }
                else{
                    errMessage['location'] = 'Visualization - UMAP - min Dist';
                    errMessage['message'] = '"min Dist" should be set';
                    return errMessage;
                }
                //nNeighbors
                if(this.DRParams['UMAP']['n_neighbors'] != ''){
                    let num = Number(this.DRParams['UMAP']['n_neighbors']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Visualization - UMAP - nNeighbors';
                        errMessage['message'] = '"nNeighbors" should be a valid number';
                        return errMessage;
                    }
                    //判断是否为正整数
                    if(num <= 0 || !Number.isInteger(num)){
                        errMessage['location'] = 'Visualization - UMAP - nNeighbors';
                        errMessage['message'] = '"nNeighbors" should be a positive integer';
                        return errMessage;
                    }
                    Params['UMAP']['n_neighbors'] = num;
                }
                else{
                    errMessage['location'] = 'Visualization - UMAP - nNeighbors';
                    errMessage['message'] = '"nNeighbors" should be set';
                    return errMessage;
                }
                //PreDR
                Params['UMAP']['PreDR'] = this.DRParams['UMAP']['PreDR']
            }
            else if(this.drMethod == 'T-SNE'){
                Params['T-SNE'] = {}
                //perplexity
                if(this.DRParams['T-SNE']['perplexity'] != ''){
                    let num = Number(this.DRParams['T-SNE']['perplexity']);
                    //判断是否为数字
                    if(isNaN(num)){
                        errMessage['location'] = 'Visualization - T-SNE - perplexity';
                        errMessage['message'] = '"perplexity" should be a valid number';
                        return errMessage;
                    }
                    //判断是否为正数
                    if(num <= 0){
                        errMessage['location'] = 'Visualization - T-SNE - perplexity';
                        errMessage['message'] = '"perplexity" should be a positive number';
                        return errMessage;
                    }
                    Params['T-SNE']['perplexity'] = num;
                }
                else{
                    errMessage['location'] = 'Visualization - T-SNE - perplexity';
                    errMessage['message'] = '"perplexity" should be set';
                    return errMessage;
                }
                //PreDR
                Params['T-SNE']['PreDR'] = this.DRParams['T-SNE']['PreDR']
            }
            return Params;
        },
    }
};
</script>

<style scoped lang="less">
.form-item {
    margin:0px;
}

</style>
