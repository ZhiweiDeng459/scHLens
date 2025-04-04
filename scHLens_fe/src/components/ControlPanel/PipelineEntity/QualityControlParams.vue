<template>
    <div>
        <!--filter Cells-->
        <el-card body-style="padding:10px" style="margin:10px 0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center">
                <b>Filter Outlying Cells</b>
                <el-switch v-model="opActive.filterCells"></el-switch>
            </div>
            <el-form :disabled="!opActive.filterCells" label-width="100px" :label-position="'left'" :model="filterCells">
                <el-form-item class="form-item" label="min Genes">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="filterCells.minGenes" size="mini" style="width: 110px;"></el-input>
                        <el-tooltip content="Minimum number of cells expressed required for a gene to pass filtering" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
                <el-form-item class="form-item" label="max Genes">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="filterCells.maxGenes" size="mini" style="width: 110px;"></el-input>
                        <el-tooltip content="Maxmum number of cells expressed required for a gene to pass filtering" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>
        
        <!--filter Genes-->
        <el-card body-style="padding:10px" style="margin:10px 0px" header-style="margin:0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center">
                <b>Filter Outlying Genes</b>
                <el-switch v-model="opActive.filterGenes"></el-switch>
            </div>
            <el-form :disabled="!opActive.filterGenes" label-width="100px" :label-position="'left'" :model="filterGenes">
                <el-form-item class="form-item" label="min Cells">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="filterGenes.minCells" size="mini" style="width: 110px;"></el-input>
                        <el-tooltip content="Minimum number of cells expressed required for a gene to pass filtering" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
                <el-form-item class="form-item" label="max Cells">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="filterGenes.maxCells" size="mini" style="width: 110px;"></el-input>
                        <el-tooltip content="Maximum number of cells expressed required for a gene to pass filtering" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
            </el-form>
        </el-card>

        <!--qc Metrics-->
        <el-card body-style="padding:10px" style="margin:10px 0px" header-style="margin:0px">
            <div slot="header" style="display:flex;justify-content:space-between;align-items:center">
                <b>MT-Gene Control</b>
                <el-switch v-model="opActive.qcMetrics"></el-switch>
            </div>
            <el-form :disabled="!opActive.qcMetrics" label-width="100px" :label-position="'left'" :model="qcMetrics">
                <el-form-item class="form-item" label="Organism">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-select v-model="qcMetrics.type" size="mini" style="width:120px">
                            <el-option 
                                v-for="o_type in organism" 
                                :key="o_type"
                                :label="o_type"
                                :value="o_type">
                            </el-option>
                        </el-select>

                        <el-tooltip content="The organism of this dataset" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item>
                <!-- <el-form-item class="form-item" label="Gene Counts">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="qcMetrics.geneCounts" size="mini" style="width: 110px;"></el-input>
                        <el-tooltip content="The number of genes with at least 1 count in a cell" placement="right">
                            <i class="el-icon-question"></i>
                        </el-tooltip>
                    </div>
                </el-form-item> -->
                <el-form-item class="form-item" label="pct Counts">
                    <div style="display:flex;justify-content:space-between;align-items:center">
                        <el-input v-model="qcMetrics.pctCounts" size="mini" style="width: 110px;"></el-input>
                        <b>%</b>
                        <el-tooltip content="Maximum proportion of total counts for a cell which are mitochondrial to pass filtering." placement="right">
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
import { Form, FormItem, Input, Tooltip, Card, Switch} from "element-ui";
import {getPipelineParamsErrorT} from "@/utils/objectTemplate";

Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);
Vue.component(Input.name, Input);
Vue.component(Tooltip.name, Tooltip);
Vue.component(Card.name, Card);
Vue.component(Switch.name, Switch);

export default {
    name: "QualityControlParams",
    props:['ExistParams'],
    data() {
        return {
            opActive:{
                //判断该操作是否启用
                filterCells:true,
                filterGenes:true,
                qcMetrics:false,
            },
            filterCells:{
                minGenes: '200',
                maxGenes: '2500',
            },
            filterGenes:{
                minCells: '3',
                maxCells: '',
            },
            qcMetrics:{
                type:'Human',
                //geneCounts: '2500',
                pctCounts: '5',
            },
            organism:[
                'Human',
                'Mouse',
            ]
        };
    },
    methods:{
        getParams(){
            /**
             * 注意要把数字字符串转为数字
             */
            let Params = {};
            let errMessage = getPipelineParamsErrorT()
            if(this.opActive.filterCells){
                Params['filterCells'] = {}

                //minGenes
                if(this.filterCells.minGenes !== ''){
                    let num = Number(this.filterCells.minGenes);
                    if(isNaN(num)){//num不合法
                        errMessage['location'] = 'Quality Control - Filter Outlying Cells - min Genes';
                        errMessage['message'] = '"min Genes" should be a valid number';
                        return errMessage;                    
                    }
                    Params['filterCells']['min_genes'] = num;
                }
                //maxGenes
                if(this.filterCells.maxGenes !== ''){
                    let num = Number(this.filterCells.maxGenes);
                    if(isNaN(num)){//num不合法
                        errMessage['location'] = 'Quality Control - Filter Outlying Cells - max Genes';
                        errMessage['message'] = '"max Genes" should be a valid number';
                        return errMessage;                    
                    }
                    Params['filterCells']['max_genes'] = Number(this.filterCells.maxGenes);
                }
            }
            if(this.opActive.filterGenes){
                Params['filterGenes'] = {}

                //minCells
                if(this.filterGenes.minCells !== ''){
                    let num = Number(this.filterGenes.minCells);
                    if(isNaN(num)){//num不合法
                        errMessage['location'] = 'Quality Control - Filter Outlying Genes - min Cells';
                        errMessage['message'] = '"min Cells" should be a valid number';
                        return errMessage;                    
                    }
                    Params['filterGenes']['min_cells'] = num;
                }
                //maxCells
                if(this.filterGenes.maxCells !== ''){
                    let num = Number(this.filterGenes.maxCells);
                    if(isNaN(num)){//num不合法
                        errMessage['location'] = 'Quality Control - Filter Outlying Genes - max Cells';
                        errMessage['message'] = '"max Cells" should be a valid number';
                        return errMessage;                    
                    }
                    Params['filterGenes']['max_cells'] = Number(this.filterGenes.maxCells);
                }
            }
            if(this.opActive.qcMetrics){
                Params['qcMetrics'] = {}
                Params['qcMetrics']['type'] = this.qcMetrics.type
                // if(this.qcMetrics.geneCounts !== '')
                //     Params['qcMetrics']['geneCounts'] = Number(this.qcMetrics.geneCounts);
                if(this.qcMetrics.pctCounts !== ''){
                    let num = Number(this.qcMetrics.pctCounts);
                    if(isNaN(num)){//num不合法
                        errMessage['location'] = 'Quality Control - MT-Gene Control - pct Counts';
                        errMessage['message'] = '"pct Counts" should be a valid number';
                        return errMessage;                    
                    }
                    if(num < 0 || num > 100){//num超出百分数范围
                        errMessage['location'] = 'Quality Control - MT-Gene Control - pct Counts';
                        errMessage['message'] = '"pct Counts" should be between 0% and 100%';
                        return errMessage;                    
                    }
                    Params['qcMetrics']['pctCounts'] = Number(this.qcMetrics.pctCounts);
                }
                else{
                    errMessage['location'] = 'Quality Control - MT-Gene Control - pct Counts';
                    errMessage['message'] = '"pct Counts" should be set';
                    return errMessage;                    
                }

            }
            return Params;
        },

    },
    watch:{
        // 'filterCells.minGenes':{
        //     handler(newValue,oldValue){
        //         if(newValue!='' && newValue !== undefined && newValue !== null){
        //             this.filterCells.maxGenes = ''
        //         }
        //     }
        // },

        // 'filterCells.maxGenes':{
        //     handler(newValue,oldValue){
        //         if(newValue!='' && newValue !== undefined && newValue !== null){
        //             this.filterCells.minGenes = ''
        //         }
        //     }
        // },

        // 'filterGenes.minCells':{
        //     handler(newValue,oldValue){
        //         if(newValue!='' && newValue !== undefined && newValue !== null){
        //             this.filterGenes.maxCells = ''
        //         }
        //     }
        // },

        // 'filterGenes.maxCells':{
        //     handler(newValue,oldValue){
        //         if(newValue!='' && newValue !== undefined && newValue !== null){
        //             this.filterGenes.minCells = ''
        //         }
        //     }            
        // }
    }
};
</script>

<style scoped lang="less">


.form-item{
    margin:0px;
}
</style>
