<template>
  <div class="AnnoRecom-container">
    <div class="AnnoRecom-title">
        <b style="font-size:18px;">{{cluster_id}} - Recommendation</b>
    </div>
    <div class="AnnoRecom-content">

        <div class="AnnoRecom-controller">
            <div style="display:flex;flex-direction:column">
                <div style="margin-bottom:10px;">
                    <b style="margin-right:12px;">Method:</b>
                    <el-select style="width:100px;margin-left:10px;" size="mini" v-model="curAnnoMethod" @change="handleMethodSelectChange">
                        <el-option
                            v-for="method in annoMethods"
                            :key="method"
                            :label="annoMethodsName[method]"
                            :value="method">
                        </el-option>
                    </el-select>
                </div>

                <div style="display:flex;">

                    <div style="margin-right:10px;">
                        <b>Gene Set:</b>
                        <el-cascader
                            size="mini"
                            v-model="curGeneSet"
                            :options="gene_sets_info"
                            @change="handleGeneSetSelectChange"
                            placeholder="Please Select A Gene Set"
                            style="margin-left:15px;width:400px;">
                        </el-cascader>
                        <el-button
                            type="primary"
                            icon="el-icon-plus"
                            size="medium"
                            style="width:26px;height:26px;padding:0px;margin-left:10px;"
                            @click="openGeneSetUploadPanel"></el-button>
                    </div>
                </div>
            </div>
            <div style="display:flex;flex-direction:column">
                <el-popover
                    width="300"
                    placement="right"
                    style="margin-bottom:10px;"
                    >
                    <div class="AnnoRecom-filter-content">
                        <!--Enricher Filter-->
                        <div class="AnnoRecom-enricher-filter-content" v-show="curAnnoMethod=='Enricher'">

                            <div style="display: flex;align-items: center;"><!--small title-->
                                <div style="border: 1px solid gray;flex:0 0 5px;height: 0px;margin-right: 6px;"></div>
                                <b style="font-size:16px">Genes Filter</b>
                                <div style="border: 1px solid gray;flex:1 1 0;height: 0px;margin-left: 6px;"></div>
                            </div>
                            <el-row type="flex" style="margin-bottom:10px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="enricher['use_logfoldchanges_threshold']"></el-checkbox>
                                </el-col>
                                <el-col :span="11">
                                    <b style="color:#666666">logfoldchanges</b>
                                    <el-tooltip style="margin-left:5px;" content='Filter by logfoldchanges. If the marker method does not calculate log fold changes, then this filter item will not be applied.' placement="right">
                                        <i class="el-icon-question"></i>
                                    </el-tooltip>
                                </el-col>
                                <el-col :span="3">&gt;</el-col>
                                <el-col :span="8"><el-input v-model="enricher['logfoldchanges_threshold']" size="mini"></el-input></el-col>
                            </el-row>
                            <!-- <el-row type="flex" style="margin-bottom:10px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="enricher['use_value_threshold']"></el-checkbox>
                                </el-col>
                                <el-col :span="11">
                                    <b style="color:#666666">value</b>
                                    <el-tooltip style="margin-left:5px;" content='Filter by expression value.' placement="right">
                                        <i class="el-icon-question"></i>
                                    </el-tooltip>
                                </el-col>
                                <el-col :span="3">&gt;</el-col>
                                <el-col :span="8"><el-input v-model="enricher['value_threshold']" size="mini" style="width:70px;"></el-input>%</el-col>
                            </el-row> -->
                            <div style="display: flex;align-items: center;margin-top: 15px;"><!--small title-->
                                <div style="border: 1px solid gray;flex:0 0 5px;height: 0px;margin-right: 6px;"></div>
                                <b style="font-size:16px">Results Filter</b>
                                <div style="border: 1px solid gray;flex:1 1 0;height: 0px;margin-left: 6px;"></div>
                            </div>
                            <el-row type="flex" style="margin-bottom:10px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="enricher['use_p_threshold']"></el-checkbox>
                                </el-col>
                                <el-col :span="11"><b style="color:#666666">FDR p-value</b></el-col>
                                <el-col :span="3">&lt;</el-col>
                                <el-col :span="8"><el-input v-model="enricher['p_threshold']" size="mini"></el-input></el-col>
                            </el-row>
                            <el-row type="flex" style="margin-bottom:0px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="enricher['use_top']"></el-checkbox>
                                </el-col>
                                <el-col :span="11">
                                    <b style="color:#666666">Number</b>
                                    <el-tooltip style="margin-left:5px;" content="Obtain the top 'Numbers' cell types with the highest combined scores." placement="right">
                                        <i class="el-icon-question"></i>
                                    </el-tooltip>
                                </el-col>
                                <el-col :span="3">&lt;</el-col>
                                <el-col :span="8"><el-input v-model="enricher['top']" size="mini"></el-input></el-col>
                            </el-row>
                        </div>

                        <!--Gsea Filter-->
                        <div class="AnnoRecom-gsea-filter-content" v-show="curAnnoMethod=='Gsea'">
                            <div style="display: flex;align-items: center;"><!--small title-->
                                <div style="border: 1px solid gray;flex:0 0 5px;height: 0px;margin-right: 6px;"></div>
                                <b style="font-size:16px">Genes Filter</b>
                                <div style="border: 1px solid gray;flex:1 1 0;height: 0px;margin-left: 6px;"></div>
                            </div>
                            <el-row type="flex" style="margin-bottom:10px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="gsea['use_logfoldchanges_threshold']"></el-checkbox>
                                </el-col>
                                <el-col :span="11">
                                    <b style="color:#666666">logfoldchanges</b>
                                    <el-tooltip style="margin-left:5px;" content='Filter by logfoldchanges. If the marker method does not calculate log fold changes, then this filter item will not be applied.' placement="right">
                                        <i class="el-icon-question"></i>
                                    </el-tooltip>
                                </el-col>
                                <el-col :span="3">&gt;</el-col>
                                <el-col :span="8"><el-input v-model="gsea['logfoldchanges_threshold']" size="mini"></el-input></el-col>
                            </el-row>
                            <div style="display: flex;align-items: center;margin-top: 15px;"><!--small title-->
                                <div style="border: 1px solid gray;flex:0 0 5px;height: 0px;margin-right: 6px;"></div>
                                <b style="font-size:16px">Results Filter</b>
                                <div style="border: 1px solid gray;flex:1 1 0;height: 0px;margin-left: 6px;"></div>
                            </div>
                            <el-row type="flex" style="margin-bottom:10px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="gsea['use_q_threshold']"></el-checkbox>
                                </el-col>
                                <el-col :span="11"><b style="color:#666666">FDR q-value</b></el-col>
                                <el-col :span="3" >&lt;</el-col>
                                <el-col :span="8"><el-input v-model="gsea['q_threshold']" size="mini"></el-input></el-col>
                            </el-row>
                            <el-row type="flex" style="margin-bottom:10px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="gsea['use_p_threshold']"></el-checkbox>
                                </el-col>
                                <el-col :span="11"><b style="color:#666666">FWER p-value</b></el-col>
                                <el-col :span="3">&lt;</el-col>
                                <el-col :span="8"><el-input v-model="gsea['p_threshold']" size="mini"></el-input></el-col>
                            </el-row>
                            <el-row type="flex" style="margin-bottom:0px;" align="middle">
                                <el-col :span="2">
                                    <el-checkbox v-model="gsea['use_top']"></el-checkbox>
                                </el-col>
                                <el-col :span="11">
                                    <b style="color:#666666">Number</b>
                                    <el-tooltip style="margin-left:5px;" content="Obtain the top 'Numbers' cell types with the highest combined scores." placement="right">
                                        <i class="el-icon-question"></i>
                                    </el-tooltip>
                                </el-col>
                                <el-col :span="3">&lt;</el-col>
                                <el-col :span="8"><el-input v-model="gsea['top']" size="mini"></el-input></el-col>
                            </el-row>
                            <!-- <el-row type="flex" style="margin-bottom:0px;" align="middle">
                                <el-col style="display:flex;justify-content:center" :span="12">
                                    <el-button style="width:80px;" size="mini" type="primary" @click="handleSetGseaFilter">Set</el-button>
                                </el-col>
                                <el-col style="display:flex;justify-content:center" :span="12">
                                    <el-button style="width:100px;" size="mini" type="danger" @click="handleCancelGseaFilter">Cancel</el-button>
                                </el-col>
                            </el-row> -->
                        </div>

                    </div>
                    <el-button slot="reference" type="primary" size="mini" style="width:100px;">Filter</el-button>
                </el-popover>
                <el-button type="success" size="mini" style="width:100px;" @click="handleRunButtonClick">Run</el-button>
            </div>
        </div>
        
        <el-tabs v-model="tabActiveName" @tab-click=handleTableClick> <!--主内容显示界面-->
            <el-tab-pane label="Table" name="Table">
                <!--enricher table -->
                <el-table
                    :data="enricher.tableData"
                    :max-height="550"
                    v-loading="enricher.loading"
                    element-loading-text="It may take some time. Please wait..."
                    @row-click="handleEnricherTableRowClick"
                    v-show="curAnnoMethod=='Enricher'"
                    border>
                    <b style="font-size:15px" slot="empty">
                        {{ emptyText_enricher }}
                    </b>
                    <el-table-column
                        prop="type"
                        label="Cell Type"
                        align="center"
                        min-width="100">
                    </el-table-column>
                    <el-table-column
                        prop="Combined Score"
                        label="Combined Score"
                        width="180"
                        align="center"
                        sortable>
                        <template slot-scope="scope">
                            <a>{{HumanFriendlyNum(scope.row['Combined Score'])}}</a>
                        </template>
                    </el-table-column>
                    <el-table-column
                        prop="Adjusted P-value"
                        label="FDR P-value"
                        align="center"
                        width="200"
                        sortable>
                        <template slot-scope="scope">
                            <a>{{HumanFriendlyNum(scope.row['Adjusted P-value'])}}</a>
                        </template>
                    </el-table-column>

                </el-table>
                <!--gsea table-->
                <el-table
                    :data="gsea.tableData"
                    :max-height="550"
                    v-loading="gsea.loading"
                    element-loading-text="It may take some time. Please wait..."
                    v-show="curAnnoMethod=='Gsea'"
                    @row-click="handleGseaTableRowClick"
                    border>
                    <b style="font-size:15px" slot="empty">
                        {{ emptyText_gsea }}
                    </b>
                    <el-table-column
                        prop="type"
                        label="Cell Type"
                        align="center"
                        min-width="100"></el-table-column>
                    <el-table-column
                        prop="NES"
                        label="NES"
                        width="90"
                        align="center"
                        sortable>
                        <template slot-scope="scope">
                            <a>{{HumanFriendlyNum(scope.row['NES'])}}</a>
                        </template>
                    </el-table-column>
                    <el-table-column
                        prop="FWER p-val"
                        label="FWER p-value"
                        align="center"
                        width="150"
                        sortable>
                        <template slot-scope="scope">
                            <a>{{HumanFriendlyNum(scope.row['FWER p-val'])}}</a>
                        </template>
                    </el-table-column>
                    <el-table-column
                        prop="FDR q-val"
                        label="FDR q-value"
                        align="center"
                        width="140"
                        sortable>
                        <template slot-scope="scope">
                            <a>{{HumanFriendlyNum(scope.row['FDR q-val'])}}</a>
                        </template>
                    </el-table-column>
                </el-table>
            </el-tab-pane>
            <el-tab-pane label="Plot" name="Plot">
                <div class="RecomResultPlotContainer" ref="RecomResultPlotContainer" v-loading="enricher.loading || gsea.loading">
                    <el-button style="position: absolute;right: 17px;top:7px;z-index: 9998;" size="mini" type="primary" icon="el-icon-download" @click="saveToFile()" circle></el-button>
                    <el-scrollbar ref="RecomResultPlot-scroll" class="recom-result-main-view">
                        <svg class="RecomResultSVG" ref="RecomResultSVG" style="background-color: white;"></svg>
                    </el-scrollbar>
                </div>
            </el-tab-pane>
        </el-tabs>

        <div>

        </div>
    </div>
    <!--基因集上传面板-->
    <UploadGeneSetsPanel ref="UploadGeneSetsPanel" :updateGeneSets="updateGeneSetsInfo" :gene_sets_info="gene_sets_info"/>

  </div>


</template>

<script>
import Vue from 'vue'
import {Row,Col,Loading,Cascader,Tabs,TabPane} from 'element-ui'
import {queryGsea,queryEnricher,queryGeneSets} from '@/utils/interface.js'
import * as d3 from "d3";
import SelfContextMenu from "@/components/SelfContextMenu"

import UploadGeneSetsPanel from "@/components/UploadGeneSetsPanel";


Vue.component(Row.name,Row)
Vue.component(Col.name,Col)
Vue.component(Cascader.name,Cascader)
Vue.component(Tabs.name,Tabs)
Vue.component(TabPane.name,TabPane)
Vue.use(Loading.directive)

export default {
    name:'AnnoRecom',
    props:['cluster_id'],
    components:{
        UploadGeneSetsPanel,
    },
    computed:{
        curData(){
            return this.$store.state.curData;
        },
        JobId(){
            return this.$store.state.JobId
        },
        ViewId(){
            return this.curData.ViewId
        },
        infoPanel(){
            return this.$store.state.infoPanel
        },
    },
    data(){
        return {
            tabActiveName:'Table',//选项卡


            curAnnoMethod:'Enricher',
            annoMethods:['Enricher','Gsea'],
            annoMethodsName:{
                'Enricher':'Enrichr',
                'Gsea':'Gsea'
            },
            // filter_visible:false,            

            //gene set
            curGeneSet:[],
            gene_sets_info:[

            ],

            //plot configure
            padding:{
                top:40,
                right:30,
                bottom:50,
                left:30
            },
            unitWidth:70,//一个unit占用的宽度
            barWidth:40,//unit中，一个bar占用的宽度

            //enricher
            init_enricher:{
                'p_threshold':0.05,
                'top':10,
                'logfoldchanges_threshold':1,
                'use_p_threshold':false,
                'use_top':false,
                'use_logfoldchanges_threshold':true,
                'tableData':[],//表格数据
                'loading':false
            },
            enricher:null,

            //gsea
            init_gsea:{
                'q_threshold':0.25,
                'p_threshold':0.05,
                'top':10,
                'logfoldchanges_threshold':1,
                'use_q_threshold':false,
                'use_p_threshold':false,
                'use_top':false,
                'use_logfoldchanges_threshold':true,
                'tableData':[], //表格数据
                'loading':false
            },
            gsea:null,

            //other
            emptyText_enricher:'',
            emptyText_gsea:'',
            
        }
    },
    methods:{
        /**
         * 外部接口
         */
        async show(curGeneSet){//在popver出现时，更新相关信息
            if(this.curGeneSet.length == 0){//如果当前没有选择基因库，那么则填入其他推荐框所选的最新基因库
                this.curGeneSet = curGeneSet;
            }

            await this.updateGeneSetsInfo();

            // this.updateRecommendResult();
        },


        /**
         * 内部函数
         */
        async updateGeneSetsInfo(){//更新基因集的相关信息
            await queryGeneSets(this.JobId).then((res)=>{
                let all_gene_sets_info = res.data;
                this.gene_sets_info = [];
                for(let org in all_gene_sets_info){
                    this.gene_sets_info.push({
                        'label':org,
                        'value':org,
                        'children':all_gene_sets_info[org].map(v=>{
                            return {
                                'value':v,
                                'label':v,
                            }
                        })
                    })
                }
            }).catch((err)=>{
                console.log('err',err)
            })

        },

        async updateRecommendResult(){//更新基因推荐的结果
            if(this.curAnnoMethod == 'Gsea'){
                await this.getGSEA();
            }
            else if(this.curAnnoMethod == 'Enricher'){
                await this.getEnricher();
            }
        },

        async getGSEA(){//获取gsea的推荐结果

            //取基因集以及进行合法性检查
            if(this.curGeneSet.length == 0){
                this.$message({
                        message:'Please select Gene Set Before Recommendation',
                        type:'error',
                        'showClose':true,
                    })
                return;
            }
            let organism = this.curGeneSet[0]
            let gene_set_name = this.curGeneSet[1]

            //取过滤参数以及进行合法性检查
            let logfoldchanges_threshold = 'all'
            let q_threshold = 'all'
            let p_threshold = 'all'
            let top = 'all'
            if(this.gsea.use_logfoldchanges_threshold){//p_threshold
                let temp = parseFloat(this.gsea.logfoldchanges_threshold)
                if(isNaN(temp)){
                    this.$message({
                        message:'Please input correct logfoldchanges threshold(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    logfoldchanges_threshold = temp
                }
            }  
            if(this.gsea.use_p_threshold){//p_threshold
                let temp = parseFloat(this.gsea.p_threshold)
                if(isNaN(temp)){
                    this.$message({
                        message:'Please input correct FWER p-value threshold(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    p_threshold = temp
                }
            }
            if(this.gsea.use_q_threshold){////q_threshold
                let temp = parseFloat(this.gsea.q_threshold)
                if(isNaN(temp)){
                    this.$message({
                        message:'Please input correct FDR q-value threshold(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    q_threshold = temp
                }
            }
            if(this.gsea.use_top){//top
                let temp = parseInt(this.gsea.top)
                if(isNaN(temp) || temp <= 0){
                    this.$message({
                        message:'Please input correct Number(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    top = temp
                }
            }

            //启动加载项
            this.gsea.loading = true;

            //发送请求
            await queryGsea(this.JobId,this.ViewId,organism,gene_set_name,this.cluster_id,logfoldchanges_threshold,q_threshold,p_threshold,top)
                .then((res)=>{
                    this.gsea.tableData = res.data.map(row=>{
                        return {
                            'type':row['Term'],
                            'NES':row['NES'],
                            'FWER p-val':(row['FWER p-val'] < Number.MIN_VALUE) ? Number.MIN_VALUE : row['FWER p-val'],
                            'FDR q-val':(row['FDR q-val'] < Number.MIN_VALUE) ? Number.MIN_VALUE : row['FDR q-val'],
                        }
                    })
                    this.gsea.loading = false
                    if (this.gsea.tableData.length > 0){
                        this.emptyText_gsea = ''
                    }
                    else{
                        this.emptyText_gsea = 'No appropriate recommended results.'
                    }
                })
                .catch((err)=>{
                    console.log('err:',err)
                    this.gsea.loading = false
                    this.gsea.tableData = []
                    this.emptyText_gsea = 'No appropriate recommended results.'

                })

        },

        async getEnricher(){//获取enricher的推荐结果
            //取基因集以及进行合法性检查
            if(this.curGeneSet.length == 0){
                this.$message({
                        message:'Please select Gene Set Before Recommendation',
                        type:'error',
                        'showClose':true,
                    })
                return;
            }
            let organism = this.curGeneSet[0]
            let gene_set_name = this.curGeneSet[1]

            //取过滤参数以及进行合法性检查
            let logfoldchanges_threshold = 'all'
            let p_threshold = 'all'
            let top = 'all'
            if(this.enricher.use_logfoldchanges_threshold){
                let temp = parseFloat(this.enricher.logfoldchanges_threshold)
                if(isNaN(temp)){
                    this.$message({
                        message:'Please input correct FDR p-value threshold(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    logfoldchanges_threshold = temp
                }
            }
            if(this.enricher.use_p_threshold){
                let temp = parseFloat(this.enricher.p_threshold)
                if(isNaN(temp)){
                    this.$message({
                        message:'Please input correct FDR p-value threshold(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    p_threshold = temp
                }
            }
            if(this.enricher.use_top){
                let temp = parseInt(this.enricher.top)
                if(isNaN(temp) || temp <= 0){
                    this.$message({
                        message:'Please input correct Number(In Filter)',
                        type:'error',
                        'showClose':true,
                    })
                    return;
                }
                else{
                    top = temp
                }
            }

            //启动加载项
            this.enricher.loading = true;

            //发送请求
            await queryEnricher(this.JobId,this.ViewId,organism,gene_set_name,this.cluster_id,logfoldchanges_threshold,p_threshold,top)
                .then((res)=>{
                    this.enricher.tableData = res.data.map(row=>{
                        return {
                            'type':row['Term'],
                            'Combined Score':row['Combined Score'],
                            'Adjusted P-value': row['Adjusted P-value'] < Number.MIN_VALUE ? Number.MIN_VALUE : row['Adjusted P-value'],
                        }
                    })
                    this.enricher.loading = false
                    if (this.enricher.tableData.length > 0){
                        this.emptyText_enricher = ''
                    }
                    else{
                        this.emptyText_enricher = 'No appropriate recommended results.'
                    }
                })
                .catch((err)=>{
                    console.log('err:',err)
                    this.enricher.loading = false
                    this.enricher.tableData = []
                    this.emptyText_enricher = 'No appropriate recommended results.'
                })

        },
        //打开基因集上传面板
        openGeneSetUploadPanel(){
            this.$refs['UploadGeneSetsPanel'].open();
        },

        /**
         * 回调函数
         */

        handleEnricherTableRowClick(row,column,event){//点击enricher表格的某一行
            this.$emit('RecommendationChosen',{
                'cluster_id':this.cluster_id,
                'type':row.type
            })

        },
        handleGseaTableRowClick(row,column,event){//点击gsea表格的某一行
            this.$emit('RecommendationChosen',{
                'cluster_id':this.cluster_id,
                'type':row.type
            })
        },
        handlePlotTypeTextClick(type){
            this.$emit('RecommendationChosen',{
                'cluster_id':this.cluster_id,
                'type':type
            })
        },
        async handleRunButtonClick(){
            await this.updateRecommendResult();
            if(this.tabActiveName == 'Plot'){
                this.reDraw()
            }
        },
        handleGeneSetSelectChange(){//gene set选择器变化
            //清空细胞类型的推荐表格
            this.gsea.tableData = []
            this.enricher.tableData = []
            //清空plot
            this.reDraw()

        },
        handleTableClick(tab){//当tab切换时触发
            setTimeout(()=>{  
                this.reDraw()
            }, 100); // 1000毫秒等于1秒  
  
        },
        handleMethodSelectChange(){//当annoMethod的select切换时出发
            this.reDraw()
        },
        /**
         * 工具函数
         */
         HumanFriendlyNum(val){//人性化显示数字 num -> str
            //过小的数字
            if(val < 0.00001)
                return val.toExponential(3) //科学计数法，保留三位有效数字
            //小于等于0的数字，但不是过小
            if(val <= 1)
                return val.toFixed(7)//保留七位小数
            //过大的数字
            if(val > 10000)
                return val.toExponential(3) //科学计数法，保留三位有效数字
            //大于0的数字，但不是过大，保留三位小数
            if(val > 0)
                return val.toFixed(3)
            //其他的，返回原本值
            return String(val)
         },

         /**
          * 画图函数
          */
         reDraw(){
            this.drawPlot(this.curAnnoMethod)
         },

         drawPlot(method){
            const self = this;
            const svg = d3.select(self.$refs['RecomResultSVG'])

            //清空视图
            svg.selectAll('*').remove()
            
            //获取根节点
            const root = svg.append('g')
            

            //获取和修饰数据
            let data = []
            if(method == 'Enricher'){
                data = self.enricher.tableData.map(v=>{
                    return {
                        'type':v['type'],
                        'value':-Math.log10(v['Adjusted P-value']),
                        'raw_value':v['Adjusted P-value']
                    }
                })
            }
            else if(method == 'Gsea'){
                data = self.gsea.tableData.map(v=>{
                    return {
                        'type':v['type'],
                        'value':-Math.log10(v['FDR q-val']),
                        'raw_value':v['FDR q-val']
                    }
                })            
            }

            //排除不合法请求
            if(data.length == 0)
                return;

            //计算宽度
            const width = self.unitWidth * data.length + self.padding.left + self.padding.right
            const height = self.$refs.RecomResultPlotContainer.clientHeight;
            svg.attr("width",width);
            svg.attr("height",height);

            //重置滚动轴
            this.$refs['RecomResultPlot-scroll'].update();

            //计算scale
            let minValue = Math.min(Math.min(...data.map(v=>v.value)),-Math.log10(1))
            let maxValue = Math.max(...data.map(v=>v.value))
            const scaleX = d3.scaleBand()
                             .domain(data.map(v=>v.type))
                             .range([self.padding.left,width - self.padding.right])
            const scaleY = d3.scaleLinear()
                             .domain([maxValue,minValue])
                             .range([self.padding.top,height - self.padding.bottom])
            
            //绘制axis
            let leftAxis = d3.axisLeft(scaleY)
            let bottomAxis = d3.axisBottom(scaleX).tickFormat(function(d){
                                    let tick_name = d
                                    // if(tick_name.length > 9){//省略名
                                    //     tick_name = tick_name.substring(0,9) + '...'
                                    // }
                                    return tick_name
                                })
            const leftAxis_g = root.append('g').call(leftAxis)
            const bottomAxis_g = root.append('g').call(bottomAxis)
            leftAxis_g.attr('transform',`translate(${self.padding.left},${0})`)
            bottomAxis_g.attr('transform',`translate(${0},${height - self.padding.bottom})`)
            
            //调整文本大小和倾斜度，并且绑定事件
            bottomAxis_g.selectAll('text')
                .attr('font-size','10px')
                .style("transform","rotate(30deg)")
                .style("transform",function(){
                    return `translate(${this.getBoundingClientRect().width * 0.5-15}px,${this.getBoundingClientRect().height * 0.5}px)` + d3.select(this).style('transform')
                })
                .on('click',(e,d)=>{
                    self.handlePlotTypeTextClick(d)
                })
                .style('cursor','pointer')
                .on("mouseover",function(e,d){ //在悬浮时突出显示
                    d3.select(this)
                        .style("fill","#3b80ee")
                })
                .on("mouseout",function(){
                    d3.select(this)
                        .style("fill",null)                        
                })

            //draw bar
            const plot = root.append('g')
            plot.selectAll('*')
                .data(data)
                .join('rect')
                .attr('x',d=>(scaleX(d.type)+0.5*self.unitWidth-0.5*self.barWidth))
                .attr('y',d=>scaleY(d.value))
                .attr('width',self.barWidth)
                .attr('height',d=>(scaleY(minValue)-scaleY(d.value)))
                .attr('fill','#0173b2')
                .on('mouseover',function(e,d){

                    let messageData = {}
                    if(method == 'Enricher'){
                        messageData['FDR P-value']=d['raw_value']
                        messageData['-lg FDR P-value']=d['value']
                    }
                    else if(method == 'Gsea'){
                        messageData['FDR q-value']=d['raw_value']
                        messageData['-lg FDR q-value']=d['value']
                    }

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData(messageData)
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)
                })
                .on('mousemove',function(e,d){

                    let messageData = {}
                    if(method == 'Enricher'){
                        messageData['FDR P-value']=d['raw_value']
                        messageData['-lg FDR P-value']=d['value']
                    }
                    else if(method == 'Gsea'){
                        messageData['FDR q-value']=d['raw_value']
                        messageData['-lg FDR q-value']=d['value']
                    }

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData(messageData)
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mouseout',function(){
                    //hidden info
                    self.infoPanel.hidden()
                })

            
            //绘制标度
            const info = root.append('g')
            info.append('text')
                .text(()=>{
                    if(method == 'Enricher'){
                        return '-lg(FDR P-value)'
                    }
                    else if(method == 'Gsea'){
                        return '-lg(FDR q-value)'
                    }
                    return ''
                })
                .attr('font-size','12px')
                .attr('x',5)
                .attr('y',20)
                .style('color','#606266')
                .style('font-family','san-serif')
            
            //重设尺寸
            svg.attr('height',root.node().getBoundingClientRect().height + self.padding.bottom)
            svg.attr('width',root.node().getBoundingClientRect().width + self.padding.right)

         },

         saveToFile(){//将svg保存为图片
            //svg保存
            const svgDOM = this.$refs['RecomResultSVG'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "bar.svg";
            a.click();
            URL.revokeObjectURL(url)

         }


    },
    watch:{
        curGeneSet(){ //选择了新的基因集，这一情况向外发送
            this.$emit('curGeneSetChange',this.curGeneSet)
        },
        curData(){ //当前数据发生变化
            //清理旧数据
            this.gsea = JSON.parse(JSON.stringify(this.init_gsea))
            this.enricher = JSON.parse(JSON.stringify(this.init_enricher))

        },
    },
    created(){
        this.gsea = JSON.parse(JSON.stringify(this.init_gsea))
        this.enricher = JSON.parse(JSON.stringify(this.init_enricher))
    }


}
</script>

<style scoped lang="less">
.AnnoRecom-container{
    display: flex;
    flex-direction: column;
    align-items: stretch;
    .AnnoRecom-title{
        padding: 5px;
        border-bottom: 3px solid gray;
        color: black;
    }
    .AnnoRecom-content{
        flex: 1 1 0;
        display: flex;
        flex-direction: column;
        .AnnoRecom-controller{
            display: flex;
            align-items: center;
            justify-content: space-between;
            margin:10px 5px;
            .AnnoRecom-filter-content{
                display: flex;
                flex-direction: column;
                .AnnoRecom-enricher-filter-content{
                    display: flex;
                    flex-direction: column;
                    padding:10px;

                }
                .AnnoRecom-gsea-filter-content{
                    display: flex;
                    flex-direction: column;
                    padding:10px;

                }

            }
        }
        .RecomResultPlotContainer{
            width: 100%;
            height: 400px;
            display: flex;
            flex-direction: column;
            align-items: stretch;
            position: relative;

            .recom-result-main-view {
                flex: 1 1;
                //overflow-x: auto;
                //overflow-y: auto;
                border-right: 2px solid rgb(200, 200, 200);
                min-height: 0;
                .RecomResultSVG{
                    
                }
                .el-scrollbar__wrap {
                    overflow: hidden;
                    width: 100%;
                }
                .is-horizontal{
                    height: 10px;
                    .el-scrollbar__thumb{
                        background-color:rgb(150, 150, 150);
                    }
                }
                .is-vertical{
                    width: 10px;
                        .el-scrollbar__thumb{
                        background-color:rgb(150, 150, 150);
                    }
            }

        }
    }
    }

    /deep/ .el-table__row{
        cursor: pointer;
    }
}
</style>