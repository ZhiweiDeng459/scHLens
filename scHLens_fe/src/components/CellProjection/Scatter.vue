<template>
    <div class="scatter-container">
        <el-scrollbar ref="scatter-scroll" class="scatter-scroller">
                <svg class="scatter" ref="scatter" xmlns="http://www.w3.org/2000/svg" style="background-color: white"></svg>
        </el-scrollbar>


        <!--右键菜单-->
        <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />

        <!--保存子数据集时，输入数据集名称的对话框-->
        <LocalDatasetPanel 
            ref="LocalDatasetPanel" 
            />

        <!--在不同投影上节点的对应视图-->
        <CorrExplorer
            ref="CorrExplorer"/>

    </div>
    
</template>

<script>
import * as d3 from "d3";
import SelfContextMenu from "@/components/SelfContextMenu"
import LocalDatasetPanel from "@/components/LocalDatasetPanel"
import CorrExplorer from "@/components/CorrExplorer"
import lasso from "@/lib/d3-lasso";
import { Popover, Table, TableColumn, Button, Loading } from 'element-ui'
import {saveSvgAsPng} from 'save-svg-png-ext'
import Vue from "vue"
import eventBus from '@/utils/eventBus.js'
import {updateGroupName,updateGroupColor,updateDeleteCells,restoreViewProjections,restoreViewLabels,fetchViewData,mergeDuplicateLabels} from '@/utils/interface'



Vue.use(Popover);
Vue.use(Table);
Vue.use(TableColumn);
Vue.use(Button);

export default {
    name: "Scatter",

    components: {
        SelfContextMenu,
        LocalDatasetPanel,
        CorrExplorer,
    },

    data() {
        return {
            // chooseColor: "#000000",
            // showColorPanel: false,
            colorIndex: 0,
            curScoreType: "Silhouette Score",
            showGlobalScore:"",
            showLocalScores:[{
                name:'test',
                score:0.0000
            }],
            menuItems:[//右键菜单项
                {
                    'name':'Select All',
                    'icon':'icons/all_select.svg',
                    'callback':()=>{
                        d3.select('.scatter-plot').selectAll('.scatter-element')
                                .classed("unchosen", false)
                                .classed("chosen", true)
                    let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map(v=>v.id)
                    this.$store.commit(`updateChosenData`, newChosenNodes);                    
                }},
                {
                    'name':'Clear Selection',
                    'callback':()=>{
                        d3.select('.scatter-plot').selectAll('.scatter-element')
                                .classed("unchosen", true)
                                .classed("chosen", false)
                        let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map(v=>v.id)
                        this.$store.commit(`updateChosenData`, newChosenNodes);
                    },
                    'icon':'icons/clear.svg'
                },
                {
                    'name':'Inverse Selection',
                    'icon':'icons/inverse_choose.svg',
                    'callback':()=>{
                        const chosen = d3.select('.scatter-plot').selectAll('.scatter-element.chosen');
                        const unchosn = d3.select('.scatter-plot').selectAll('.scatter-element.unchosen')
                        chosen.classed("unchosen", true)
                            .classed("chosen", false) 
                        unchosn.classed("unchosen", false)
                            .classed("chosen", true)  
                        let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map(v=>v.id)
                        this.$store.commit(`updateChosenData`, newChosenNodes);
                    },
                },
                // {
                //     'name':'Hide Chosen Elements',
                //     'icon':'icons/hide.svg',
                //     'callback':()=>{
                //         const chosen = d3.select('.scatter-plot').selectAll('.scatter-element.chosen');
                //         //清空被选择点
                //         chosen.classed("unchosen", true)
                //             .classed("chosen", false) 
                //         let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map(v=>v.id)
                //         this.$store.commit(`updateChosenData`, newChosenNodes);
                //         chosen.classed('hidden',true)


                //     }
                // },
                // {
                //     'name':'Show Hidden Elements',
                //     'icon':'icons/show.svg',
                //     'callback':()=>{
                //         const hidden = d3.selectAll('.hidden')
                //         hidden.classed('hidden',false)
                //     }
                // },
                {//合并同名点
                    'name':'Merge Duplicate Labels',
                    'icon':'icons/mergeLabels.svg',
                    'callback':()=>{
                        const loading = Loading.service({ fullscreen: true });

                        mergeDuplicateLabels(this.JobId,this.ViewId)
                            .then((response)=>{
                                fetchViewData(this.JobId,this.ViewId)
                                    .then((fetch_response)=>{
                                        let newViewData = fetch_response.data
                                        console.log('newViewData',newViewData)
                                        this.$message({
                                            'message':'The merge process finished',
                                            'type':'success',
                                            'showClose':true,
                                        })
                                        //提交新的数据
                                        this.$store.commit("updateViewData",newViewData)
                                        //切换到更新的视图
                                        this.$store.commit("toggleCurData", this.ViewId);

                                        loading.close();

                                    })
                                    .catch((fetch_err)=>{
                                        console.log(fetch_err)
                                        this.$message({
                                            'message':'A error ocurred in running process of the merge process.',
                                            'type':'error',
                                            'showClose':true,
                                        })  
                                        loading.close();
                                    })    
                            })
                            .catch((err)=>{
                                console.log(err)
                                this.$message({
                                    'message':'A error ocurred in running process of the merge process.',
                                    'type':'error',
                                    'showClose':true,
                                })  
                                loading.close();
                            })
                    }
                }, 
                {//恢复标签
                    'name':'Restore the Labels',
                    'icon':'icons/restore.svg',
                    'callback':()=>{
                        const loading = Loading.service({ fullscreen: true ,text:'It may take several minutes, please wait...'});

                        restoreViewLabels(this.JobId,this.ViewId)
                            .then((response)=>{
                                fetchViewData(this.JobId,this.ViewId)
                                    .then((fetch_response)=>{
                                        let newViewData = fetch_response.data
                                        console.log('newViewData',newViewData)
                                        this.$message({
                                            'message':'The restore process finished',
                                            'type':'success',
                                            'showClose':true,
                                        })
                                        //提交新的数据
                                        this.$store.commit("updateViewData",newViewData)
                                        //切换到更新的视图
                                        this.$store.commit("toggleCurData", this.ViewId);

                                        loading.close();

                                    })
                                    .catch((fetch_err)=>{
                                        console.log(fetch_err)
                                        this.$message({
                                            'message':'A error ocurred in running process of the restore process.',
                                            'type':'error',
                                            'showClose':true,
                                        })  
                                        loading.close();
                                    })


                            })
                            .catch((err)=>{
                                console.log(err)
                                this.$message({
                                    'message':'A error ocurred in running process of the restore process.',
                                    'type':'error',
                                    'showClose':true,
                                })  
                                loading.close();
                            })
                    }
                },
                {//恢复投影
                    'name':'Restore the Projection',
                    'icon':'icons/restore.svg',
                    'callback':()=>{
                        const loading = Loading.service({ fullscreen: true ,text:'It may take several minutes, please wait...'});

                        restoreViewProjections(this.JobId,this.ViewId)
                            .then((response)=>{
                                fetchViewData(this.JobId,this.ViewId)
                                    .then((fetch_response)=>{
                                        let newViewData = fetch_response.data
                                        console.log('newViewData',newViewData)
                                        this.$message({
                                            'message':'The restore process finished',
                                            'type':'success',
                                            'showClose':true,
                                        })
                                        //提交新的数据
                                        this.$store.commit("updateViewData",newViewData)
                                        //切换到更新的视图
                                        this.$store.commit("toggleCurData", this.ViewId);

                                        loading.close();

                                    })
                                    .catch((fetch_err)=>{
                                        console.log(fetch_err)
                                        this.$message({
                                            'message':'A error ocurred in running process of the restore process.',
                                            'type':'error',
                                            'showClose':true,
                                        })  
                                        loading.close();
                                    })


                            })
                            .catch((err)=>{
                                console.log(err)
                                this.$message({
                                    'message':'A error ocurred in running process of the restore process.',
                                    'type':'error',
                                    'showClose':true,
                                })  
                                loading.close();
                            })
                    }
                },
                {//查询对应点
                    'name':'Discover Correspondence',
                    'icon':'icons/find.svg',
                    'callback':()=>{
                        this.$refs['CorrExplorer'].openExplorer()
                    }
                }, 
                {//删除被选择点
                    'name':'Delete Chosen Elements',
                    'icon':'icons/delete.svg',
                    'callback':()=>{
                        const loading = Loading.service({ fullscreen: true });

                        updateDeleteCells(this.JobId,this.ViewId,this.chosenData)
                            .then((response)=>{
                                /**
                                 * 返回的是更新数据
                                 * 
                                 * 更新要点：
                                 * 
                                 * 1.清空chosenCells，因为已经被删除
                                 * 2.迭代式更新以下curData中的以下数据
                                 *  a. cellData
                                 *  b. globalScores和localScores
                                 *  c. groups
                                 *  d. MK
                                 *  e. raw_embedding_range
                                 * 3. 重绘当前视图 
                                 * 
                                 * 注意考虑细胞删除完的情况
                                 * 
                                 */

                                let self = this

                                //清空选择数据
                                d3.select('.scatter-plot').selectAll('.scatter-element')
                                        .classed("unchosen", true)
                                        .classed("chosen", false)
                                let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map(v=>v.id)
                                console.log('newChosenNodes:',newChosenNodes)
                                this.$store.commit(`updateChosenData`, newChosenNodes);

                                let updateTree = response.data;
                                
                                function updateByTree(cur){
                                    let findedData = self.dataList.find((v)=>{
                                        return v.ViewId == cur.ViewId
                                    })
                                    findedData.MK = cur.updateData.MK;
                                    findedData.cellData = cur.updateData.cellData;
                                    findedData.groups = cur.updateData.groups;
                                    findedData.localScores = cur.updateData.localScores;
                                    findedData.globalScores = cur.updateData.globalScores;
                                    for(let child of cur.children){
                                        updateByTree(child)
                                    }
                                }
                                updateByTree(updateTree)

                                loading.close()
                                //
                            })
                            .catch((err)=>{
                                console.log(err)
                                loading.close()
                            })
                    }
                },
                {
                    'name':'Create Local Plot',
                    'icon':'icons/drill.svg',
                    'callback':()=>{
                        eventBus.$emit('callPipelineConfig','local');
                    }
                },
                {
                    'name':'Export the Chosen Cells',
                    'icon':'icons/save.svg',
                    'callback':()=>{
                        this.$refs['LocalDatasetPanel'].openDialog();//打开对话框
                    }
                },
                {
                    'name':'Save this Image',
                    'icon':'icons/save_as_image.svg',
                    'callback':()=>{
                        this.saveToFile();
                    }
                }

            ]

        };
    },

    computed: {
        curData(){
            return this.$store.state.curData;
        },
        dataList(){
            return this.$store.state.dataList
        },
        groups(){
            return this.curData.groups;
        },
        cellData(){
            return this.curData.cellData;
        },
        globalScores() {
            return this.curData.globalScores;
        },
        localScores(){
            return this.curData.localScores;
        },
        repaintTag(){
            return this.$store.state.repaintTag;
        },
        chosenData(){
            return this.$store.state.curData.chosenData;
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
    watch: {
        cellData() {
            if (this.cellData === undefined || this.cellData === null) return;
                this.reDraw()
        },
        groups: {
            //监视
            deep: true,
            handler() {
                if (this.groups === "undefined" || this.groups === "null") return;
                this.reDraw()
            },
        },
        'repaintTag.Scatter':{
            handler(){
                //TODO更新分数
                this.reDraw();
            }
        },


        chosenData(newValue,oldValue){
            //设置散点央视
            d3.select('.scatter-plot').selectAll('.scatter-element')
                .classed('chosen',(d,i)=>{
                    if(newValue.find(v=>v == d.id) !== undefined)
                        return true
                    else return false
                })
                .classed('unchosen',d=>{
                    if(newValue.find(v=>v == d.id) !== undefined)
                        return false
                    else return true  
                })
            //修改chosen text
            d3.select(this.$refs['scatter']).select('.chosen-num-text').text(`Chosen Cells : ${newValue.length}`)
        }
        
    },
    methods: {
        updateScore(){
            //global

            //local
            this.showLocalScores.length = 0;
            for(let group of this.groups){
                this.showLocalScores.push({'name':group['name'],'score':this.localScores[group.id]})
            }
        },

        drawScatter() {

            const self = this;

            const pointArr = JSON.parse(JSON.stringify(this.cellData));
            const labels = JSON.parse(JSON.stringify(this.groups));


            const svg = d3.select(`.scatter`);
            svg.selectAll("*").remove();
            const legend = svg.append('g')
            legend.classed('scatter-legend',true);
            const scatter = svg.append('g')
            scatter.classed('scatter-plot',true);



            /**
             * 计算各个尺寸信息
             */
             const height = this.$refs.scatter.clientHeight;   //svg总高度

            const scatterHeight = height;
            const scatterWidth = height * 0.9
            const scatterPadding = 60;

            // const legendWidth = width * 0.1;
            const legendHeight = height;
            const legendPadding = 20;


            /**
             * draw legend
             */

            legend.attr("transform", `translate(${scatterWidth},${0})`)

            //generate legend unit circle
            legend
                .selectAll(`.scatter-legend-item`)
                .data(labels)
                .enter()
                .append("circle")
                .classed("scatter-legend-item", true)
                .attr("cx", 0)
                .attr("cy", 0)
                .attr("fill", (d) => d.color)
                .attr("r", 5)
                .style("cursor", "pointer")
                .style('stroke-width',3)
                .on("click", function (event,d) {//选择/取消选择一类点
                    console.log('d:',d)
                    if(this.classList.contains('chosen')){//选中 -> 不选中
                        this.classList.add('unchosen')
                        this.classList.remove('chosen')                      
                        const unchosen =  d3.select('.scatter-plot').selectAll(`.scatter-element.label${d.id.replaceAll(" ","，").replaceAll(",","，").replaceAll("+","，").replaceAll("~","，").replaceAll(">","，").replaceAll(".","，")}`)
                            .classed("unchosen", true)
                            .classed("chosen", false)
                    }
                    else{//不选中 -> 选中
                        this.classList.add('chosen')
                        this.classList.remove('unchosen') 
                        const chosen = d3.select('.scatter-plot').selectAll(`.scatter-element.label${d.id.replaceAll(" ","，").replaceAll(",","，").replaceAll("+","，").replaceAll("~","，").replaceAll(">","，").replaceAll(".","，")}`)
                            .classed("unchosen", false)
                            .classed("chosen", true)
                    }
                    //更新数据
                    let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map(v=>v.id)
                    self.$store.commit(`updateChosenData`, newChosenNodes);

                })
            
            //generate legend unit text
            legend
                .selectAll(`.scatter-legend-item-text`)
                .data(labels)
                .enter()
                .append("text")
                .attr("x", 0)
                .attr("y", 0)
                .text((d) => d.name)
                .attr('dominant-baseline','central')
                .style('font-family','YaHei')
                .style('font-weight','bold')
                .classed("scatter-legend-item-text", true)
                .style("cursor", "pointer")
                .on("mouseover",function(){ //在悬浮时突出显示
                    d3.select(this)
                        .style("fill","#3b80ee")
                })
                .on("mouseout",function(){
                    d3.select(this)
                        .style("fill","black")                        
                })
                .on('click',function(e,d){
                    svg.selectAll('.scatter-legend-item').filter(v=>v.id==d.id).dispatch('click')
                })

            //motify legend position
            let legendUnitHeight = 30;
            let legendNumInAColumn = Math.floor(1.0 * (legendHeight - 2 * legendPadding) / legendUnitHeight) //一列的legend unit数目
            let columnNum = Math.ceil(1.0 * labels.length / legendNumInAColumn) //有多少列
            let maxTextWidth = new Array(columnNum).fill(0); //数组保存了每列的最长文本长度
            legend.selectAll('.scatter-legend-item-text')//calculate max text width
                  .each(function(d,i){
                    let curColumnIndex = Math.floor(1.0 * i / legendNumInAColumn)
                    if(this.getBoundingClientRect().width > maxTextWidth[curColumnIndex]){
                        maxTextWidth[curColumnIndex] = this.getBoundingClientRect().width;
                    }
                  })
            let legendUnitWidth = maxTextWidth.map(v=>v + 20)
            let legendAccumulateWitdh = [0] //legend每列的累计宽度，第一项为0，最后一项为所有unit的宽度和，因此legend的列数多一
            for(let i = 0;i < legendUnitWidth.length;i++){
                legendAccumulateWitdh.push(legendUnitWidth[i] + legendAccumulateWitdh[i])
            }
            let legendWidth = legendAccumulateWitdh[legendAccumulateWitdh.length-1] + 2 * legendPadding

            legend.selectAll('.scatter-legend-item')
                  .attr('cx',(d,i)=>{
                    let curColumnIndex = Math.floor(1.0 * i / legendNumInAColumn)
                    return legendPadding + legendAccumulateWitdh[curColumnIndex] + 10;
                  })
                  .attr('cy',(d,i)=>{
                    return legendPadding + (i % legendNumInAColumn) * legendUnitHeight + 0.5 * legendUnitHeight;
                  })
            legend.selectAll('.scatter-legend-item-text')
                  .attr('x',(d,i)=>{
                    let curColumnIndex = Math.floor(1.0 * i / legendNumInAColumn)
                    return legendPadding + legendAccumulateWitdh[curColumnIndex] + 20;
                  })
                  .attr('y',(d,i)=>{
                    return legendPadding + (i % legendNumInAColumn) * legendUnitHeight + 0.5 * legendUnitHeight;
                  })
            



            /**
             * draw scatter
             */

            let minX = this.curData.raw_embedding_range['x'][0]
            let maxX = this.curData.raw_embedding_range['x'][1]
            let minY = this.curData.raw_embedding_range['y'][0]
            let maxY = this.curData.raw_embedding_range['y'][1]

            const posXScale = d3
                .scaleLinear()
                .domain([minX, maxX])
                .range([scatterPadding, scatterWidth -  scatterPadding]);
            const posYScale = d3
                .scaleLinear()
                .domain([minY, maxY])
                .range([scatterPadding, scatterHeight - scatterPadding]);


            const circles = scatter
                .selectAll("circle")
                .data(pointArr)
                .enter()
                .append("circle")
                .attr("cx", (d) => posXScale(d.pos[0]))
                .attr("cy", (d) => posYScale(d.pos[1]))
                .attr("r", Math.min(2000.0 / pointArr.length,3))
                .attr("fill", (d) => this.groups.find(e=>d.group == e.id).color)
                .attr("fill-opacity",0.7)
                .attr("stroke",(d) => this.groups.find(e=>d.group == e.id).color)
                .attr("stroke-width","1px")
                .attr("id", (d) => d.id)
                .attr("class", (d) => `label${d.group.replaceAll(" ","，").replaceAll(",","，").replaceAll("+","，").replaceAll("~","，").replaceAll(">","，").replaceAll(".","，")}`) //将空格去掉
                .classed('scatter-element',true)
                .classed("unchosen", function(d,i){
                    if(self.chosenData.find(v=>v == d.id) !== undefined)
                        return false;
                    return true;
                })
                .classed("chosen",function(d,i){
                    if(self.chosenData.find(v=>v == d.id) !== undefined)
                        return true;
                    return false;
                })
                .on('label_select',(e)=>{
                    let label = e.detail.label
                })
                .on('mouseover',function(e,d){

                    let anno = self.groups.find(g=>{
                        return g.id == d.group
                    }).name

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'id':d.id,'group':anno})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mousemove',function(e,d){
                    let anno = self.groups.find(g=>{
                        return g.id == d.group
                    }).name

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'id':d.id,'group':anno})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mouseout',function(){
                    //hidden info
                    self.infoPanel.hidden()

                })



            //add score
            const score = svg.append("text")
                            .attr("x",10)
                            .attr("y",20)
                            .text(this.showGlobalScore)
                            .style('font-family','YaHei')
                            .style('font-weight','bold')

            
            // add chosen num text
            svg.append("text")
                .classed('chosen-num-text',true)
                .attr("x",10)
                .attr("y",height - 10)
                .text(`Chosen Cells : ${self.chosenData.length}`)
                .style('font-family','YaHei')
                .style('font-weight','bold')
                .style('font-size',13)



            //set canvas size
            let parentWidth = d3.select('.scatter-container').node().clientWidth
            let width = scatterWidth + legendWidth > parentWidth ? (scatterWidth + legendWidth) : parentWidth
            svg.attr('width',width)


            //lasso
            var lasso_start = () => {
            };
            var lasso_draw = () => {
                ls.possibleItems().classed("unchosen", false).classed("chosen", true)
            };
            var lasso_end = () => {
                ls.selectedItems().classed("unchosen", false).classed("chosen", true).classed("possible", false).dispatch("chosen");
                
                let newChosenNodes = d3.select('.scatter-plot').selectAll('.scatter-element.chosen').data().map((v)=>{
                    return v.id;
                })
                self.$store.commit(`updateChosenData`, newChosenNodes);
            };
            var ls = lasso()
                .closePathSelect(true)
                .closePathDistance(100)
                .items(circles)
                .targetArea(svg)
                .on("start", lasso_start)
                .on("draw", lasso_draw)
                .on("end", lasso_end);
            ls.closePathDistance(2000);
            svg.call(ls);

        },



        menuMounted(_this,root,parent) {
        
        },

        saveToFile(){


            /**
             * 保存前
             */
             d3.select(this.$refs['scatter']).select('.chosen-num-text').style('display','none') //隐藏选择的细胞数



            /**
             * 保存视图为文件
             */
            // saveSvgAsPng(this.$refs.scatter, "scatter.png");
            const svgDOM = this.$refs['scatter'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "scatter.svg";
            a.click();
            URL.revokeObjectURL(url)

            /**
             * 结尾处理
             */
             d3.select(this.$refs['scatter']).select('.chosen-num-text').style('display',null) //显示选择的细胞数


        },
        reDraw(){
            eventBus.$emit("CellProjectionViewRefreshingStart")

            if(this.curData === undefined || this.curData === null)
                return ;
            if(this.curData.cellData === undefined || this.curData.cellData === null){
                return ;
            }


            if(this.curData.cellData.length == 0){//处理数据量为0的情况
                d3.select('.scatter').selectAll('*').remove();//删除所有投影图图像
            }
            else{
                //加载分数
                if(this.globalScores[this.curScoreType] != 'No score')
                    this.showGlobalScore = this.curScoreType + ': ' + this.globalScores[this.curScoreType].toFixed(3);
                else
                    this.showGlobalScore = this.curScoreType + ': '  + this.globalScores[this.curScoreType]
                this.updateScore();
                //绘图
                this.drawScatter();
            }
            eventBus.$emit("CellProjectionViewRefreshingClose")

            this.$refs['scatter-scroll'].update()//更新滚动轴长度
        }
    },
    mounted(){
        this.reDraw()
    }
};
</script>

<style lang="less" scoped>
.scatter-container {
    width: 100%;
    height: 100%;
    display: inline-block;
    display: flex;
    position: relative;
    align-items: stretch;

    .scatter-scroller{
        flex: 1 1;
        position:relative; 
    }

    .color-picker{
        position: absolute;
        z-index: 10;
        // top:5% !important;
        // left:100% !important;
    }
    .group-info{
        right:10%;
        top:1%;
        position:absolute;
        z-index: 1000;
    }

    .scatter {
        height: 100%;
    }
    .scatter-legend-item-text {
        font-size: 12px;
    }
}


circle {
    .chosen,
    .possible {
        fill-opacity: 1;
        z-index: 100;
        stroke: black;

    }
    .unchosen {
        fill-opacity: 0.9;
    }

}

.hidden {
    display: none;
}

.lasso {
    path {
        fill-opacity: 0.6;
        stroke: rgb(64, 169, 255);
        stroke-width: 2px;
    }
    .lasso {
        .drawn {
            fill-opacity: 0.05;
        }
        .loop_close {
            fill: none;
            stroke-dasharray: 4, 4;
        }
        .lasso .origin {
            fill: #3399ff;
            fill-opacity: 0.5;
        }
    }
}

/**
* scroll
*/

/deep/ .el-scrollbar__wrap{
    overflow: hidden;
    position: relative;
    width:100%;
    .el-scrollbar__view{
        height: 100%;
    }
}
/deep/ .is-horizontal{
    height: 10px!important;
    .el-scrollbar__thumb{
        background-color:rgb(150, 150, 150);
    }
}
/deep/ .is-vertical{
    width: 10px;
    display: none;
    .el-scrollbar__thumb{
        background-color:rgb(150, 150, 150);
    }
}

</style>
