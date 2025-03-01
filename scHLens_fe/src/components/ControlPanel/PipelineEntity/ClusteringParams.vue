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
import { Form, FormItem, Input, Select, Option, Radio, Tooltip } from "element-ui";

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);

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
            if(this.clusterMethod == 'leiden'){
                Params['leiden'] = {}
                if(this.clusterParams['leiden']['resolution'] != ''){
                    Params['leiden']['resolution'] = Number(this.clusterParams['leiden']['resolution']);
                }
                if(this.clusterParams['leiden']['resolution'] != ''){
                    Params['leiden']['n_neighbors'] = Number(this.clusterParams['leiden']['n_neighbors']);
                }
            }
            else if(this.clusterMethod == 'kmeans'){
                Params['kmeans'] = {}
                if(this.clusterParams['kmeans']['n_clusters'] != ''){
                    Params['kmeans']['n_clusters'] = Number(this.clusterParams['kmeans']['n_clusters']);
                }
            }
            // else if(this.clusterMethod == 'louvain'){
            //     Params['louvain'] = {}
            //     if(this.clusterParams['louvain']['resolution'] != ''){
            //         Params['louvain']['resolution'] = Number(this.clusterParams['louvain']['resolution']);
            //     }
            //     if(this.clusterParams['louvain']['resolution'] != ''){
            //         Params['louvain']['n_neighbors'] = Number(this.clusterParams['louvain']['n_neighbors']);
            //     }
            // }
            else if(this.clusterMethod == 'sc3s'){
                Params['sc3s'] = {}
                if(this.clusterParams['sc3s']['n_clusters'] != ''){
                    Params['sc3s']['n_clusters'] = Number(this.clusterParams['sc3s']['n_clusters']);
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
