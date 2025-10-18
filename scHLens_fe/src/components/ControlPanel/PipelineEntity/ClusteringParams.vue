<template>
    <div>

        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="clusterMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in clusterOptions" :key="item.value" :label="item.label"
                        :value="item.value"></el-option>
                </el-select>
            </div>

            <!--k-means-->
            <div v-if="clusterMethod == 'kmeans'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['kmeans']">
                    <el-form-item class="form-item" label="MultiK">
                        <div style="display:flex;justify-content:space-between;align-items:center;height:40px">
                            <el-switch v-model="clusterParams['kmeans']['auto_number']"></el-switch>
                            <el-tooltip
                                content="Automated Cluster Number Recommendation using MultiK. Enabling this option will override 'n_clusters'."
                                placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="n_clusters">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['kmeans']['n_clusters']" size="mini" style="width: 110px;"
                                :disabled="kmeans_nclusters_disabled"></el-input>
                            <el-tooltip content="The target number of clusters" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <div style="padding:5px 5px">
                        <i class="el-icon-warning" style="color:orange"></i>
                        <a style="word-break: normal; overflow-wrap: normal;">
                            Since MultiK performs multiple rounds of clustering and consistency calculations, it can be
                            quite time-consuming. Therefore, for large datasets, we recommend identifying an appropriate
                            number of clusters through multiple explorations with visualization views to save computation
                            time.
                        </a>
                    </div>
                </el-form>
            </div>


            <!--leiden-->
            <div v-if="clusterMethod == 'leiden'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['leiden']">
                    <el-form-item class="form-item" label="MultiK">
                        <div style="display:flex;justify-content:space-between;align-items:center;height:40px">
                            <el-switch :disabled="true"></el-switch>
                            <el-tooltip
                                content="Automated Cluster Number Recommendation using MultiK. Enabling this option will override 'n_clusters'."
                                placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <div style="padding:5px 5px">
                        <i class="el-icon-warning" style="color:orange"></i>
                        <a style="word-break: normal; overflow-wrap: normal;">
                            Since MultiK only provides recommendations for the number of clusters, it cannot be used to
                            suggest the resolution parameter for the Leiden algorithm. Therefore, the automatic
                            recommendation feature of MultiK is not available for Leiden. To use this function, please
                            switch the clustering method to k-means or sc3s.
                        </a>
                    </div>
                    <el-form-item class="form-item" label="nNeighbors">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['leiden']['n_neighbors']" size="mini"
                                style="width: 110px;"></el-input>
                            <el-tooltip placement="right">
                                <div slot="content">
                                    <div style="display: flex;align-items: center;">
                                        <span>The size of local neighborhood</span>
                                        <el-button
                                            @click="tutorial_jump('/tutorial/startaAnalysisPipeline/2/t_Leiden_nNeighbors', 't_Leiden_nNeighbors')"
                                            size="mini" type="primary"
                                            style="margin-left: 10px;height: 20px;width:90px;padding: 0px;display: flex;align-items: center;justify-content: center;">
                                            More Details
                                        </el-button>
                                    </div>
                                </div>
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="resolution">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['leiden'].resolution" size="mini"
                                style="width: 110px;"></el-input>
                            <el-tooltip placement="right">
                                <div slot="content">
                                    <div style="display: flex;align-items: center;">
                                        <span>resolution</span>
                                        <el-button
                                            @click="tutorial_jump('/tutorial/startaAnalysisPipeline/2/t_Leiden_resolution', 't_Leiden_resolution')"
                                            size="mini" type="primary"
                                            style="margin-left: 10px;height: 20px;width:90px;padding: 0px;display: flex;align-items: center;justify-content: center;">
                                            More Details
                                        </el-button>
                                    </div>
                                </div>
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
            <div v-if="clusterMethod == 'sc3s'">
                <el-form label-width="100px" :label-position="'left'" :model="clusterParams['sc3s']">
                    <el-form-item class="form-item" label="MultiK">
                        <div style="display:flex;justify-content:space-between;align-items:center;height:40px">
                            <el-switch v-model="clusterParams['sc3s']['auto_number']"></el-switch>
                            <el-tooltip
                                content="Automated Cluster Number Recommendation using MultiK. Enabling this option will override 'n_clusters'."
                                placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <el-form-item class="form-item" label="n_clusters">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-input v-model="clusterParams['sc3s'].n_clusters" size="mini" style="width: 110px;"
                                :disabled="sc3s_nclusters_disabled"></el-input>
                            <el-tooltip content="The target number of clusters" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                    <div style="padding:5px 5px">
                        <i class="el-icon-warning" style="color:orange"></i>
                        <a style="word-break: normal; overflow-wrap: normal;">
                            Since MultiK performs multiple rounds of clustering and consistency calculations, it can be
                            quite time-consuming. Therefore, for large datasets, we recommend identifying an appropriate
                            number of clusters through multiple explorations with visualization views to save computation
                            time.
                        </a>
                    </div>
                </el-form>
            </div>

        </el-card>
    </div>
</template>

<script>
import Vue from "vue";
import { getPipelineParamsErrorT } from "@/utils/objectTemplate";
import { Form, FormItem, Input, Select, Option, Radio, Tooltip, MessageBox } from "element-ui";
import eventBus from "@/utils/eventBus.js"

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
            kmeans_nclusters_disabled: false, //控制kmeans的n_clusters输入框是否禁用
            sc3s_nclusters_disabled: false, //控制sc3s的n_clusters输入框是否禁用
            clusterParams: {
                'kmeans': {
                    'n_clusters': '8',
                    'auto_number': false,
                },
                'leiden': {
                    'resolution': '1',
                    'n_neighbors': '10',
                },
                // 'louvain':{
                //     'resolution':'1',
                //     'n_neighbors':'10',
                // },
                'sc3s': {
                    'n_clusters': '8',
                    'auto_number': false,
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
    methods: {
        getParams() {
            /**
             * 注意要把数字字符串转为数字
             */
            let Params = {};
            let errMessage = getPipelineParamsErrorT();

            if (this.clusterMethod == 'leiden') {

                Params['leiden'] = {}

                //resolution
                if (this.clusterParams['leiden']['resolution'] != '') {
                    let num = Number(this.clusterParams['leiden']['resolution']);
                    if (!Number.isNaN(num)) {
                        if (num > 0) {
                            Params['leiden']['resolution'] = num
                        }
                        else {//如果resolution小于等于0，则报错
                            errMessage['location'] = 'Clustering - leiden - resolution';
                            errMessage['message'] = '"resolution" should be a positive number';
                            return errMessage;
                        }

                    }
                    else {//如果resolution不是数字，则报错
                        errMessage['location'] = 'Clustering - leiden - resolution';
                        errMessage['message'] = '"resolution" should be a valid positive number';
                        return errMessage;
                    }
                }
                else {//如果没有设置resolution，则报错
                    errMessage['location'] = 'Clustering - leiden - resolution';
                    errMessage['message'] = '"resolution" should be set';
                    return errMessage;
                }

                //nNeighbors
                if (this.clusterParams['leiden']['n_neighbors'] != '') {
                    let num = Number(this.clusterParams['leiden']['n_neighbors']);
                    if (!Number.isNaN(num)) {
                        if (num > 0 && Number.isInteger(num)) {
                            Params['leiden']['n_neighbors'] = num
                        }
                        else {//如果n_neighbors小于等于0或者不是整数，则报错
                            errMessage['location'] = 'Clustering - leiden - nNeighbors';
                            errMessage['message'] = '"nNeighbors" should be a positive integer';
                            return errMessage;
                        }
                    }
                    else {//如果n_neighbors不是数字，则报错
                        errMessage['location'] = 'Clustering - leiden - nNeighbors';
                        errMessage['message'] = '"nNeighbors" should be a valid positive integer';
                        return errMessage;
                    }
                }
                else {//如果没有设置n_neighbors，则报错
                    errMessage['location'] = 'Clustering - leiden - nNeighbors';
                    errMessage['message'] = '"nNeighbors" should be set';
                    return errMessage;
                }

            }

            else if (this.clusterMethod == 'kmeans') {
                Params['kmeans'] = {}
                //n_clusters
                if (this.clusterParams['kmeans']['n_clusters'] != '') {
                    let num = Number(this.clusterParams['kmeans']['n_clusters']);
                    if (!Number.isNaN(num)) {
                        if (num > 0 && Number.isInteger(num)) {
                            Params['kmeans']['n_clusters'] = num
                        }
                        else {//如果n_clusters小于等于0或者不是整数，则报错
                            errMessage['location'] = 'Clustering - k-means - n_clusters';
                            errMessage['message'] = '"n_clusters" should be a positive integer';
                            return errMessage;
                        }
                    }
                    else {//如果n_clusters不是数字，则报错
                        errMessage['location'] = 'Clustering - k-means - n_clusters';
                        errMessage['message'] = '"n_clusters" should be a valid positive integer';
                        return errMessage;
                    }
                }
                else {//如果没有设置n_clusters，则报错
                    errMessage['location'] = 'Clustering - k-means - n_clusters';
                    errMessage['message'] = '"n_clusters" should be set';
                    return errMessage;
                }
                Params['kmeans']['auto_number'] = this.clusterParams['kmeans']['auto_number']
            }

            else if (this.clusterMethod == 'sc3s') {
                Params['sc3s'] = {}
                if (this.clusterParams['sc3s']['n_clusters'] != '') {
                    let num = Number(this.clusterParams['sc3s']['n_clusters']);
                    if (!Number.isNaN(num)) {
                        if (num > 0 && Number.isInteger(num)) {
                            Params['sc3s']['n_clusters'] = num
                        }
                        else {//如果n_clusters小于等于0或者不是整数，则报错
                            errMessage['location'] = 'Clustering - sc3s - n_clusters';
                            errMessage['message'] = '"n_clusters" should be a positive integer';
                            return errMessage;
                        }
                    }
                    else {//如果n_clusters不是数字，则报错
                        errMessage['location'] = 'Clustering - sc3s - n_clusters';
                        errMessage['message'] = '"n_clusters" should be a valid positive integer';
                        return errMessage;
                    }
                }
                else {//如果没有设置n_clusters，则报错
                    errMessage['location'] = 'Clustering - sc3s - n_clusters';
                    errMessage['message'] = '"n_clusters" should be set';
                    return errMessage;
                }
                Params['sc3s']['auto_number'] = this.clusterParams['sc3s']['auto_number']
            }
            return Params;
        },
        tutorial_jump(index, id) {
            eventBus.$emit('TutorialJump', index, id)
        }
    },
    watch: {
        'clusterParams.kmeans.auto_number': {
            handler(newVal) {
                if (newVal) {
                    //启动自动推荐，则禁用n_clusters输入框
                    this.kmeans_nclusters_disabled = true
                }
                else {
                    this.kmeans_nclusters_disabled = false
                }
            }
        },
        'clusterParams.sc3s.auto_number': {
            handler(newVal) {
                if (newVal) {
                    //启动自动推荐，则禁用n_clusters输入框
                    this.sc3s_nclusters_disabled = true
                }
                else {
                    this.sc3s_nclusters_disabled = false
                }
            }
        },
    }
};
</script>

<style scoped lang="less">.form-item {
    margin: 0px;
}</style>
