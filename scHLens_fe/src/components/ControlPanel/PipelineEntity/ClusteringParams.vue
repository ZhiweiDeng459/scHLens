<template>
    <div>
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="clusterMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in clusterOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
            </div >

            <!--k-means-->
            <div v-if="clusterMethod=='kmeans'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['kmeans']">
                    <el-form-item class="form-item" label="n_clusters">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['kmeans']['n_clusters']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The final number of clusters" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>


            <!--leiden-->
            <div v-if="clusterMethod=='leiden'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['leiden']">
                    <el-form-item class="form-item" label="nNeighbors">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['leiden']['n_neighbors']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The size of local neighborhood" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="resolution">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['leiden'].resolution" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="resolution" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>

            <!--louvain-->
            <!-- <div v-if="clusterMethod=='louvain'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['louvain']">
                    <el-form-item class="form-item" label="nNeighbors">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['louvain']['n_neighbors']" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The size of local neighborhood" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="resolution">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['louvain'].resolution" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="resolution" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div> -->

            <!--sc3s-->
            <div v-if="clusterMethod=='sc3s'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['sc3s']">
                    <el-form-item class="form-item" label="n_clusters">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['sc3s'].n_clusters" size="mini" style="width: 110px;"></el-input>
                            <el-tooltip content="The final number of clusters" placement="right">
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
import {getPipelineParamsErrorT} from "@/utils/objectTemplate";
import { Form, FormItem, Input, Select, Option, Radio, Tooltip, MessageBox } from "element-ui";

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);
Vue.component(MessageBox.name, MessageBox);

export default {
    name: "ClusterParams",
    data() {
        return {
            clusterMethod: "leiden",
            clusterParams:{
                'kmeans':{
                    'n_clusters':'8'
                },
                'leiden':{
                    'resolution':'1',
                    'n_neighbors':'10',
                },
                // 'louvain':{
                //     'resolution':'1',
                //     'n_neighbors':'10',
                // },
                'sc3s':{
                    'n_clusters':'8'
                },
            },
            clusterOptions: [
                {
                    value: "kmeans",
                    label: "k-means",
                },
                {
                    value: "leiden",
                    label: "leiden",
                },
                // {
                //     value: "louvain",
                //     label: "louvain",
                // },
                {
                    value: "sc3s",
                    label: "sc3s",
                },
            ],
        };
    },
    methods:{
        getParams(){
            /**
             * 注意要把数字字符串转为数字
             */
            let Params = {};
            let errMessage = getPipelineParamsErrorT();

            if(this.clusterMethod == 'leiden'){
                
                Params['leiden'] = {}
                
                //resolution
                if(this.clusterParams['leiden']['resolution'] != ''){
                    let num = Number(this.clusterParams['leiden']['resolution']);
                    if(!Number.isNaN(num)){
                        if(num > 0){
                            Params['leiden']['resolution'] = num
                        }
                        else{//如果resolution小于等于0，则报错
                            errMessage['location'] = 'Clustering - leiden - resolution';
                            errMessage['message'] = '"resolution" should be a positive number';
                            return errMessage;
                        }
                        
                    }
                    else{//如果resolution不是数字，则报错
                        errMessage['location'] = 'Clustering - leiden - resolution';
                        errMessage['message'] = '"resolution" should be a valid positive number';
                        return errMessage;
                    }
                }
                else{//如果没有设置resolution，则报错
                    errMessage['location'] = 'Clustering - leiden - resolution';
                    errMessage['message'] = '"resolution" should be set';
                    return errMessage;
                }

                //nNeighbors
                if(this.clusterParams['leiden']['n_neighbors'] != ''){
                    let num = Number(this.clusterParams['leiden']['n_neighbors']);
                    if(!Number.isNaN(num)){
                        if(num > 0 && Number.isInteger(num)){
                            Params['leiden']['n_neighbors'] = num
                        }
                        else{//如果n_neighbors小于等于0或者不是整数，则报错
                            errMessage['location'] = 'Clustering - leiden - nNeighbors';
                            errMessage['message'] = '"nNeighbors" should be a positive integer';
                            return errMessage;
                        }
                    }
                    else{//如果n_neighbors不是数字，则报错
                        errMessage['location'] = 'Clustering - leiden - nNeighbors';
                        errMessage['message'] = '"nNeighbors" should be a valid positive integer';
                        return errMessage;
                    }
                }
                else{//如果没有设置n_neighbors，则报错
                    errMessage['location'] = 'Clustering - leiden - nNeighbors';
                    errMessage['message'] = '"nNeighbors" should be set';
                    return errMessage;
                }

            }
            
            else if(this.clusterMethod == 'kmeans'){
                Params['kmeans'] = {}
                //n_clusters
                if(this.clusterParams['kmeans']['n_clusters'] != ''){
                    let num = Number(this.clusterParams['kmeans']['n_clusters']);
                    if(!Number.isNaN(num)){
                        if(num > 0 && Number.isInteger(num)){
                            Params['kmeans']['n_clusters'] = num
                        }
                        else{//如果n_clusters小于等于0或者不是整数，则报错
                            errMessage['location'] = 'Clustering - k-means - n_clusters';
                            errMessage['message'] = '"n_clusters" should be a positive integer';
                            return errMessage;
                        }
                    }
                    else{//如果n_clusters不是数字，则报错
                        errMessage['location'] = 'Clustering - k-means - n_clusters';
                        errMessage['message'] = '"n_clusters" should be a valid positive integer';
                        return errMessage;
                    }
                }
                else{//如果没有设置n_clusters，则报错
                    errMessage['location'] = 'Clustering - k-means - n_clusters';
                    errMessage['message'] = '"n_clusters" should be set';
                    return errMessage;
                }
            }

            else if(this.clusterMethod == 'sc3s'){
                Params['sc3s'] = {}
                if(this.clusterParams['sc3s']['n_clusters'] != ''){
                    let num = Number(this.clusterParams['sc3s']['n_clusters']);
                    if(!Number.isNaN(num)){
                        if(num > 0 && Number.isInteger(num)){
                            Params['sc3s']['n_clusters'] = num
                        }
                        else{//如果n_clusters小于等于0或者不是整数，则报错
                            errMessage['location'] = 'Clustering - sc3s - n_clusters';
                            errMessage['message'] = '"n_clusters" should be a positive integer';
                            return errMessage;
                        }
                    }
                    else{//如果n_clusters不是数字，则报错
                        errMessage['location'] = 'Clustering - sc3s - n_clusters';
                        errMessage['message'] = '"n_clusters" should be a valid positive integer';
                        return errMessage;
                    }
                }
                else{//如果没有设置n_clusters，则报错
                    errMessage['location'] = 'Clustering - sc3s - n_clusters';
                    errMessage['message'] = '"n_clusters" should be set';
                    return errMessage;
                }
            }
            return Params;
        }
    }
};
</script>

<style scoped lang="less">

.form-item {
    margin:0px;
}

</style>
