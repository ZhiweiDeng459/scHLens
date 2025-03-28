<template>
    <div class="Tutorial-container">
            <el-menu
                :default-active="$route.path"  
                mode="vertical"
                background-color="#24292f"
                text-color="#fff"
                active-text-color="#fff"
                class="tutorial-menu"
                @open="handleSubMenuClick"
                @close="handleSubMenuClick"
                router>
                <el-scrollbar ref="t-scroll" class="t-scroll">
                    <el-menu-item class="menu-item" index="/tutorial/intro">
                        <span class="menu-item-text-1" slot="title">Introduction</span>
                    </el-menu-item>
                    <el-submenu class="menu-item" index="/tutorial/enterjob">
                        <span class="menu-item-text-1" slot="title">Enter a Job</span>
                        <el-menu-item index="/tutorial/enterjob/1" @click="scrollTo('t_WhatisJob')">What is Job?</el-menu-item>
                        <el-menu-item index="/tutorial/enterjob/2" @click="scrollTo('t_WhatisJobID')">What is Job ID?</el-menu-item>
                        <el-menu-item index="/tutorial/enterjob/3" @click="scrollTo('t_CreateLoadJob')">Create a New Job & Load an Existing Job & Use Samples</el-menu-item>
                        <el-menu-item index="/tutorial/enterjob/4" @click="scrollTo('t_ExportImportJob')">Export the Job & Import the Job</el-menu-item>
                        <el-menu-item index="/tutorial/enterjob/5" @click="scrollTo('t_DeleteJob')">Delete Job</el-menu-item>
                    </el-submenu>
                    <el-submenu class="menu-item" index="/tutorial/selectDataset">
                        <span class="menu-item-text-1" slot="title">Select a Dataset</span>
                        <el-menu-item index="/tutorial/selectDataset/1" @click="scrollTo('t_SelectSample')">Select Sample Dataset</el-menu-item>
                        <el-menu-item index="/tutorial/selectDataset/2" @click="scrollTo('t_AddDataset')">Add your own Dataset</el-menu-item>
                    </el-submenu>
                    <el-menu-item class="menu-item" index="/tutorial/startaAnalysisPipeline">
                        <span class="menu-item-text-1" slot="title">Start an Analysis Pipeline</span>
                    </el-menu-item>
                    <el-submenu class="menu-item" index="/tutorial/analysiswithViews">
                        <span class="menu-item-text-1" slot="title">Analysis with Views</span>
                        <el-submenu class="menu-item" index="/tutorial/analysiswithViews/1">
                            <div slot="title" style="width:100%;height: 100%;" @click="scrollTo('t_Views')">
                                <span class="menu-item-text-2"  >Views</span>
                            </div>
                            
                            <el-menu-item index="/tutorial/analysiswithViews/1/1" @click="scrollTo('t_CellPrjectionView')">Cell Projection View</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/1/2" @click="scrollTo('t_GeneProjectionView')">Gene Projection View</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/1/3" @click="scrollTo('t_DiffView')">Distribution of Gene Expression View</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/1/4" @click="scrollTo('t_MarkerGeneView')">Marker Gene View</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/1/5" @click="scrollTo('t_HierarchyNavigation1')">Hierarchy Navigation</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/1/6" @click="scrollTo('t_Controller')">Controller Bar</el-menu-item>
                        </el-submenu>
                        <el-submenu class="menu-item" @click="scrollTo('t_Interactions')" index="/tutorial/analysiswithViews/2">
                            <div slot="title" style="width:100%;height: 100%;" @click="scrollTo('t_Interactions')">

                            <span class="menu-item-text-2">Interactions</span>
                            </div>
                            <el-menu-item index="/tutorial/analysiswithViews/2/1" @click="scrollTo('t_SelectCells')">Select Cells</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/2" @click="scrollTo('t_DeleteCells')">Delete Cells</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/3" @click="scrollTo('t_ExportCells')">Export Cells</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/4" @click="scrollTo('t_MergeDuplicateLabels')">Merge Duplicate Labels(Annotations)</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/5" @click="scrollTo('t_Recovery')">Recovery</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/6" @click="scrollTo('t_SwitchGenes')">Switch Current Genes</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/7" @click="scrollTo('t_EditGroup')">Edit group color & Export group markers</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/8" @click="scrollTo('tAnnotation')">Annotation</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/9" @click="scrollTo('t_SaveImages')">Save the images</el-menu-item>
                            <el-menu-item index="/tutorial/analysiswithViews/2/10" @click="scrollTo('t_Refreshing')">Refreshing Views</el-menu-item>
                        </el-submenu>
                    </el-submenu>
                    <el-submenu class="menu-item" index="/tutorial/hierarchicalExploration">
                        <span class="menu-item-text-1" slot="title">Hierarchical Exploration</span>
                        <el-menu-item index="/tutorial/hierarchicalExploration/1" @click="scrollTo('t_CrateLocalPlot')">Create Local Plot & Find Good Pattern</el-menu-item>
                        <el-menu-item index="/tutorial/hierarchicalExploration/2" @click="scrollTo('t_HierarchyNavigation2')">Hierarchy Navigation</el-menu-item>
                        <el-menu-item index="/tutorial/hierarchicalExploration/3" @click="scrollTo('t_DiscoverCorrespondence')">Discover the correspondence of cells across different nodes</el-menu-item>
                        <el-menu-item index="/tutorial/hierarchicalExploration/4" @click="scrollTo('t_Merge')">Merge the Local Results</el-menu-item>
                    </el-submenu>
                </el-scrollbar>
            </el-menu>
            <router-view/>
    </div>
</template>

<script>
import Vue from 'vue';


export default {
    name: "Tutorial",
    components: {

    },
    computed: {
        curData() {
            return this.$store.state.curData;
        },
        dataList(){
            return this.$store.state.dataList;
        },
        pipelineEntityArray(){
            return this.$store.state.pipelineEntityArray;
        },
        JobId(){
            return this.$store.state.JobId
        }
    },
    methods:{
        scrollTo(id) {
            // this.$router.push({ hash: id }).catch(err => {});
            this.$nextTick(()=>{
                let element = document.getElementById(id);
                if (element) {
                    setTimeout(() => {
                        element.scrollIntoView({ behavior: 'smooth' });
                    }, 0);
                }

                //重置滚动轴
                this.$refs['t-scroll'].update();
            })
            

        },
        handleSubMenuClick(index,indexPath){
            this.$router.push(index)


            this.$refs['t-scroll'].update();

        }

    }


};
</script>

<style scoped lang="less">
    .Tutorial-container{
        background-color:#F5F5F5;
        width: 100%;
        height: 100%;
        display: flex;
        overflow: hidden;
        .tutorial-menu{
            flex: 0 0 440px;
            overflow-x:hidden;
           .menu-item{
                margin:0px 20px 0px 0px;
                .menu-item-text-1{ //一级标题
                    font-size:18px;
                }
                .menu-item-text-2{ //二级标题
                    font-size:15px;
                }
                /deep/ .el-submenu__icon-arrow{
                    color:white;
                }
           }

        }
        /deep/ .t-scroll{
            overflow:hidden;
            min-height: 0;
            min-width: 0;
            padding:20px 0px;
            .el-scrollbar__wrap {
                    overflow: hidden;
                }
            .is-horizontal{
                    height: 10px;
                    display: none;
                    .el-scrollbar__thumb{
                        
                        background-color:rgb(150, 150, 150);
                    }
                }
            .is-vertical{
                    width: 10px;
                    display: none;
                .el-scrollbar__thumb{
                    background-color:rgb(150, 150, 150);
                }
            }

        }
    }
    
</style>