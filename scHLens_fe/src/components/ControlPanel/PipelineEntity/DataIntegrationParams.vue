<template>
    <div>
        <!-- Methods -->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Method</b>
                <el-select v-model="DIMethod" placeholder="Select..." size="mini" style="width:150px">
                    <el-option v-for="item in DIMethodOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                </el-select>
            </div >
            <!--Ingest-->
            <div v-if="DIMethod=='Ingest'">
                <el-form label-width="100px" :label-position="'left'">
                    <el-form-item class="form-item" label="Embedding">
                        <div style="display:flex;justify-content:space-between;align-items:center">
                            <el-select v-model="DIParams['Ingest']['embeddingMethod']" placeholder="Select..." size="mini" style="width:130px">
                                <el-option v-for="item in IngestEmbeddingOptions" :key="item.value" :label="item.label" :value="item.value"></el-option>
                            </el-select>
                            <el-tooltip content="Minimum Dist" placement="right">
                                <i class="el-icon-question"></i>
                            </el-tooltip>
                        </div>
                    </el-form-item>
                </el-form>
            </div>
            <!--Scanorama-->
            <div v-if="DIMethod=='Scanorama'">
                <el-form label-width="100px" :label-position="'left'">
                </el-form>
            </div>
            <!--Harmony-->
            <div v-if="DIMethod=='Harmony'">
                <el-form label-width="100px" :label-position="'left'">
                </el-form>
            </div>
        </el-card>

        <!-- Dataset -->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center;">
                <b>Dataset</b>
                <el-select v-model="dataset" placeholder="Select..." size="mini" style="width:120px">
                    <el-option v-for="item in dataSetOptions" :key="item.index" :label="item.label" :value="item" :disabled="(integrationDatasets.map(e=>e.index).indexOf(item.index) != -1) || item.index == ref_dataset.index"></el-option>
                </el-select>
                <el-button type="primary" icon="el-icon-plus" style="width:30px;height:30px;padding:0px" @click="addDataSet"></el-button>
            </div>
            <el-form>
                <el-form-item v-for="item in integrationDatasets" :key="item.index" style="margin:0px;padding:0px">
                    <div class="integration-entity">
                        {{item.label}}
                        <div>
                            <el-popover
                                placement="right"
                                width="400"
                                trigger="click">
                                <QualityControlParams ref="qcParams"/>
                                <el-button type="primary" slot="reference" circle icon="el-icon-s-tools" style="padding:3px;margin:10px 4px"></el-button>
                            </el-popover>
                            <el-button type="info" circle icon="el-icon-delete-solid" style="padding:3px;margin:10px 0px" @click="deleteDataSet(item.index)"></el-button>
                        </div>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>
    </div>

</template>

<script>
import Vue from "vue"
import { Form, FormItem, Input, Select, Option, Radio, Tooltip, Card} from "element-ui";
import QualityControlParams from "@/components/ControlPanel/PipelineEntity/QualityControlParams";
Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Select.name, Select);
Vue.component(Option.name, Option);
Vue.component(Radio.name, Radio);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Card.name, Card);
export default {
    name:"DataIntegrationParams",
    props:['ref_dataset'],
    components:{
        QualityControlParams,
    },
    data(){
        return {
          dataset:null,
          DIParams:{
            'Ingest':{
                embeddingMethod:'PCA'
            },
          },
          DIMethod:"Ingest",
          DIMethodOptions:[
            {
                value: "Ingest",
                label: "Ingest",
            },
            {
                value: "Scanorama",
                label: "Scanorama",
            },
            {
                value: "Harmony",
                label: "Harmony",
            },
          ],
          integrationDatasets:[],
          IngestEmbeddingOptions:[
            {
                value: "PCA",
                label: "PCA",
            },
            {
                value: "UMAP",
                label: "UMAP",
            },
          ],
        }
    },
    computed:{
        dataSetOptions(){
            return this.$store.state.dataSetOptions;
        }
    },
    watch:{
        ref_dataset(){
            this.integrationDatasets.splice(0,this.integrationDatasets.length)
        }
    },
    methods:{
        addDataSet(){
            if(this.dataset != null){
                this.integrationDatasets.push(this.dataset);
                this.dataset=null
            }
        },
        deleteDataSet(index){
            this.integrationDatasets.splice(this.integrationDatasets.map(e=>e.index).indexOf(index),1);
        },
        getParams(){
            /**
             * 注意要把数字字符串转为数字
             */
            let Params = {};
            Params['Datasets'] = []
            for(let i = 0;i < this.integrationDatasets.length;i++){
                Params['Datasets'].push({
                    'Dataset':this.integrationDatasets[i].value,
                    'qcParams':this.$refs['qcParams'][i].getParams()
                })
            }
            if(this.DIMethod == 'Ingest'){
                Params['Method'] = {}
                Params['Method']['Ingest'] = this.DIParams['Ingest']
            }
            else if(this.DIMethod == 'Scanorama'){
                Params['Method'] = {}
                Params['Method']['Scanorama'] = {}
            }
            else if(this.DIMethod == 'Harmony'){
                Params['Method'] = {}
                Params['Method']['Harmony'] = {}
            }
            return Params;
        }        
    },
}
</script>s

<style lang="less">
.integration-entity{
    border-bottom: 3px solid #F5F5F5;
    display: flex;
    flex-direction: row;
    justify-content: space-between;
}
</style>