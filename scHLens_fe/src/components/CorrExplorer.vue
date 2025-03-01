
<!--节点之间的对应关系的查找组件-->

<template>
    <el-dialog
        width="1400px"
        :visible.sync="showDialog"
        :show-close="true"
        :lock-scroll="true"
        :modal="true"
        :close-on-click-modal="false"
        title="Correspondence Explorer"
    >
        <div class="corr-explorer-container">

            <div ref="tree-container" class="tree-container">
                <div class="tree-title-container">
                    <b>Tree Structure</b>
                </div>
                <svg ref="tree-plot"  class="tree-plot"></svg>
                <ScatterSketch ref="scatter-sketch"/>

            </div>
            <div class="scatter-conatiner">
                <div class="scatter-title-container">
                    <b>Correspondence Overview</b>
                </div>
                <svg ref="scatter-plot" class="scatter-plot"></svg>

            </div>
        </div>
    </el-dialog>

</template>

<script>

import * as d3 from "d3";
import ScatterSketch from "@/components/ScatterSketch"
import curNodeImg from '@/assets/icons/curNode.svg'
import NodeImg from '@/assets/icons/Node.svg'


export default {

    

    name: "CorrExplorer",
    
    components:{
        ScatterSketch,
    },


    data(){
        return {
            /**
             * 关键数据
             */
            showDialog:false,//对话框是否显示
            reference:null,//参考节点

            /**
             * 配置文件
             */
            //Scatter plot
            NoScatterViewInfoTitleSize:28,
            NoScatterViewInfoContentSize:19,
            ScatterPadding:{
                'top' : 40,
                'left' : 30,
                'right' : 30,
                'bottom' : 60,
            },
            backgroundCellColor:'#217FC6', //背景细胞的颜色
            localCellColor:'#DDA620', //局部细胞的颜色
            legend_unit_height:60,//legend一个单元的高度
            legend_circle_size:10,//legend中圆的大小
            legend_text_size:15,//legend中文本的大小
            
            //Tree plot
            TreePadding:{
                'top' : 40,
                'left' : 30,
                'right' : 30,
                'bottom' : 60,
            },
            sketchWidth:150,
            sketchHeight:150,
            arrowSize:13,
            nodeTextSize:15,

        }
    },

    computed:{
        curData(){
            return this.$store.state.curData;
        },
        projectTree() {
            return this.$store.state.projectTree;
        },
        dataList(){
            return this.$store.state.dataList;
        },
        chosenData(){
            return this.$store.state.curData.chosenData;
        },
    },

    methods:{
        
        /**
         * 外部接口
         */

        openExplorer(){//打开对话框
            const self = this;
            self.showDialog = true;
            this.$nextTick(()=>{
                self.reDraw();
            })
        },

        closeExplorer(){//关闭对话框
            const self = this;
            self.showDialog = true;
        },
        
        /**
         * 内部函数
         */
        reDrawTree(){
            const svg = d3.select(this.$refs['tree-plot'])
            const self = this;
            //clear
            svg.selectAll('*').remove();
            //设定尺寸
            const width = self.$refs['tree-plot'].getBoundingClientRect().width
            const height = self.$refs['tree-plot'].getBoundingClientRect().height


            const tree = d3.tree().size([width - self.TreePadding.left - self.TreePadding.right,height - self.TreePadding.top - self.TreePadding.bottom])
            //准备数据
            const rootNode = self.projectTree['root'];
            let nodeStack = [rootNode];
            let drawData = {};
            let drawStack = [drawData] 
            while(nodeStack.length != 0){
                let curNode = nodeStack.pop();
                let curDraw = drawStack.pop();

                curDraw['name'] = curNode.id;
                
                if(curNode['children'].length != 0){  
                    curDraw['children'] = []
                    for(let i = 0; i < curNode['children'].length;i++){
                        let temp = {}
                        curDraw['children'].push(temp);
                        drawStack.push(temp);
                        nodeStack.push(curNode['children'][i]);
                    }
                }
            }

            let root = tree(d3.hierarchy(drawData));

            if(root.height == 0){//单层
                root.y = 0.5 * (height - self.TreePadding.top - self.TreePadding.bottom)
            }


            const nodeSize = Math.min(50,0.5 * (height - self.TreePadding.top - self.TreePadding.left) / (root.height + 1))
            const SelectRingRadius = 0.5 * nodeSize


            //select Ring
            const selectRing = svg.append('g')
            root.each(d => {
                selectRing
                .append('circle')
                .attr('r',SelectRingRadius)
                .attr('cx',d.x + self.TreePadding.left)
                .attr('cy',d.y + self.TreePadding.top)
                .attr('stroke','black')
                .attr('stroke-width',6)
                .attr('fill',()=>{
                    return 'white'
                })
                .classed('unselectNode',function(){
                    if(d['data']['name'] != self.reference)
                        return true;
                    else{
                        return false;
                    }
                })
            });


            //entity
            const Entity = svg.append('g')
            root.each(d=>{
                Entity.append('image')
                    .attr('width',nodeSize)
                    .attr('height',nodeSize)
                    .attr('x',function(){
                        return d.x +  self.TreePadding.left - this.getBoundingClientRect().width * 0.5
                    })
                    .attr('y',function(){
                        return d.y +  self.TreePadding.top - this.getBoundingClientRect().height * 0.5;
                    })
                    .attr('href',function(){
                        return (d['data']['name'] == self.curData['ViewId'] ?  curNodeImg : NodeImg)
                    })
                    .style('cursor','pointer')
                    .on('click',function(){//选择/取消选择元素
                        self.reference = d['data']['name']
                        self.reDraw();
                    })
                    .on("contextmenu", function (event) {//切换视图
                        event.preventDefault();
                        self.reference = d['data']['name']
                        self.reDraw();
                        event.stopPropagation();
                        
                    })
                    .on("mouseover",function(event){//显示略缩图
                        
                        const sketch = self.$refs['scatter-sketch'];

                        //显示
                        sketch.show();

                        //整理数据
                        let targetData = self.dataList.find(v=>{
                            return v['ViewId'] == d['data']['name']
                        });
                        let sketchData = targetData.cellData.map(v=>{
                            return {
                                'x':v.pos[0],
                                'y':v.pos[1],
                                'color':targetData.groups.find(e=>v.group == e.id).color
                            }
                        })

                        //设置尺寸、位置与绘图
                        sketch.setSize(self.sketchWidth,self.sketchHeight)

                        const containerWidth = self.$refs['tree-container'].getBoundingClientRect().width
                        const containerHeight = self.$refs['tree-container'].getBoundingClientRect().height
                        let pos = [event.offsetX + 20,event.offsetY + 20]
                        if(pos[0] + self.sketchWidth + 12 >= containerWidth){
                            pos[0] = event.offsetX - self.sketchWidth - 7
                        }
                        if(pos[1] + self.sketchHeight + 12 >= containerHeight){
                            pos[1] = event.offsetY - self.sketchHeight - 7
                        }
                        sketch.setPos(pos[0],pos[1])

                        sketch.draw(sketchData,targetData['raw_embedding_range'])
                    })
                    .on('mousemove',function(event){  

                        const sketch = self.$refs['scatter-sketch'];
                        const containerWidth = self.$refs['tree-container'].getBoundingClientRect().width
                        const containerHeight = self.$refs['tree-container'].getBoundingClientRect().height
                        let pos = [event.offsetX + 10,event.offsetY + 10]
                        if(pos[0] + self.sketchWidth + 12 >= containerWidth){
                            pos[0] = event.offsetX - self.sketchWidth - 7
                        }
                        if(pos[1] + self.sketchHeight + 12 >= containerHeight){
                            pos[1] = event.offsetY - self.sketchHeight - 7
                        }
                        sketch.setPos(pos[0],pos[1])
                    })
                    .on("mouseout",function(event){//隐藏略缩图
                        const sketch = self.$refs['scatter-sketch'];
                        sketch.hidden();
                    })
            })
            //修复切换时移动略缩图不消失的bug
            svg.on("mouseout",function(event){
                const sketch = self.$refs['scatter-sketch'];
                sketch.hidden();
            })



            
            //EntityText
            const EntityText = svg.append('g')
            root.each(d=>{
                    EntityText.append('text')
                        .datum(d['data']['name'])
                        .text((text)=>text)
                        .style('font-size',self.nodeTextSize)
                        .attr('x',function(){
                            return d.x + self.TreePadding.left - this.getBoundingClientRect().width * 0.5
                        })
                        .attr('y',function(){
                            return d.y + self.TreePadding.top + self.nodeTextSize * 2 + 10
                        })
                        .style('font-family','YaHei')
                        .style('font-weight','bold')
            })


            //arrow
            svg.append('defs')
               .append('marker')
               .attr('id','corr_explorer_arrow')
               .attr('markerWidth',self.arrowSize)
               .attr('markerHeight',self.arrowSize)
               .attr('refX',0)
               .attr('refY',0.5 * self.arrowSize)
               .attr("markerUnits","userSpaceOnUse")
               .attr('orient', 'auto')
               .append('path')
               .attr('fill', '#434343')
               .attr('d', `M 0,0 L ${self.arrowSize},${0.5*self.arrowSize} L 0,${self.arrowSize}`)


            //link
            const links = svg.append('g')
                .selectAll('path')
                .data(root.links())
                .join('path')
                .attr('fill','none')
                .attr('stroke','#666666')
                .attr('stroke-width',5)
                .attr('marker-end','url(#corr_explorer_arrow)')
                .attr('d',function(data){
                    const link = d3.link(d3.curveBumpY)
                                   .source(function(d){
                                        const sourceText = EntityText.selectAll('text')
                                            .filter((text)=>{
                                                return d.source.data.name == text
                                            })
                                            .node()
                                        
                                        return {
                                            'x':d.source.x,
                                            'y':parseFloat(sourceText.getAttribute('y')) + 10
                                        }
                                   })
                                   .target(function(d){
                                        
                                        return {
                                            'x':d.target.x,
                                            'y':d.target.y + self.TreePadding.top - self.arrowSize - 25,
                                        }
                                   })
                                   .x(d=>d.x + self.TreePadding.left)
                                   .y(d=>d.y)
                    return link(data)
                })




        },

        reDrawScatter(){
            const svg = d3.select(this.$refs['scatter-plot'])
            const self = this;
            //clear
            svg.selectAll('*').remove();
            //设定尺寸
            const width = self.$refs['scatter-plot'].getBoundingClientRect().width;
            const height = self.$refs['scatter-plot'].getBoundingClientRect().height;
            const scatterWidth = 0.8 * width;
            const legendwidth = 0.2 * width;

            let ref_data = undefined;
            let cur_data = self.curData;
            let idAvailableFlag = false //id是否合法
            let inhertAvailableFlag = false; //是否有继承关系


            //检查reference id是否能找到对应的数据
            if(self.reference !== undefined && self.reference !== null){//检索reference data
                ref_data = self.dataList.find((e)=>{
                    if(e['ViewId'] == self.reference){
                        return true;
                    }
                    return false;
                })  
            }
            if(ref_data === undefined || ref_data === null){
                idAvailableFlag = false
            }
            else{
                idAvailableFlag = true
            }

            //校验是否有继承关系
            function check_children(root,find_id){
                    if (root.id == find_id){
                        inhertAvailableFlag = true;
                        return;
                    }
                    else{
                        for(let child of root['children']){
                            check_children(child,find_id)
                        }
                    }
                    return;
                }

            if(ref_data !== undefined && ref_data !== null){
                let ref_node = ref_data['TreeNode']
                let cur_node = cur_data['TreeNode']
                check_children(ref_node,cur_node.id)
                // check_children(cur_node,ref_node.id)
            }

            if(!idAvailableFlag || !inhertAvailableFlag){//不能合法显示scatter
                svg.append('text')
                    .text('No Available View')
                    .style('font-size',self.NoScatterViewInfoTitleSize)
                    .attr('text-anchor','middle')
                    .attr('dominant-baseline','middle')
                    .attr('x','50%')
                    .attr('y','40%')
                    .style('font-family','YaHei')
                    .style('font-weight','bold')

                svg.append('text')
                    .text('Please select the ancestor of the current node')
                    .style('font-size',self.NoScatterViewInfoContentSize)
                    .attr('text-anchor','middle')
                    .attr('dominant-baseline','middle')
                    .attr('x','50%')
                    .attr('y','50%')
                    .style('font-family','YaHei')
                return ;
            }

            //获取数据
            let backgroundCells = null;
            let localIds = null;
            let localMap = {}
            backgroundCells = JSON.parse(JSON.stringify(ref_data.cellData))
            // localIds = cur_data.map((v)=>v.id)
            if(self.chosenData.length != 0){
                localIds = JSON.parse(JSON.stringify(self.chosenData))
            }
            else{
                localIds = cur_data.cellData.map(v=>v.id)
            }
            //映射集合
            for(let id of localIds){
                localMap[id] = true
            }
            //计算半径
            let radius = 4000.0 / backgroundCells.length
            //计算范围
            let range = {
                'x':[
                    Math.min(...backgroundCells.map(v=>v.pos[0])),
                    Math.max(...backgroundCells.map(v=>v.pos[0])),
                ],
                'y':[
                    Math.min(...backgroundCells.map(v=>v.pos[1])),
                    Math.max(...backgroundCells.map(v=>v.pos[1])),
                ]
            }
            //计算scale
            const xScale = d3.scaleLinear()
                .domain([...range.x])
                .range([radius + self.ScatterPadding.left,scatterWidth-radius-self.ScatterPadding.right])
            const yScale = d3.scaleLinear()
                .domain([...range.y])
                .range([radius + self.ScatterPadding.top,height-radius-self.ScatterPadding.bottom])
            //绘图
            const plot = svg.append('g')
                .selectAll('*')
                .data(backgroundCells)
                .join('circle')
                .attr('cx',d=>xScale(d.pos[0]))
                .attr('cy',d=>yScale(d.pos[1]))
                .attr('fill',d=>{
                    if(d.id in localMap){
                        return self.localCellColor
                    }
                    else{
                        return self.backgroundCellColor
                    }
                })
                .attr('r',radius)
            const legend = svg.append('g')
            const legend_unit = legend
                            .selectAll('*')
                            .data([
                                {'name':'Reference'},
                                {'name':self.chosenData.length == 0?'Current':'Chosen'}
                            ])
                            .join('g')
                            .attr('transform',(d,i)=>{
                                return `translate(${scatterWidth},${(i + 1) * self.legend_unit_height})`
                            })
            legend_unit.append('circle')
                        .attr('r',7)
                        .attr('fill',(d)=>{
                            if(d.name == 'Reference')
                                return self.backgroundCellColor
                            else
                                return self.localCellColor
                        })
                        .attr('cx',10)
                        .attr('cy',self.legend_unit_height * 0.5)
            
            legend_unit.append('text')
                       .text(d=>d.name)
                       .attr('dominant-baseline','middle')
                       .attr('x',10 + 7 + 5)
                       .attr('y',self.legend_unit_height * 0.5)
                       .style('font-family','YaHei')
                       .style('font-weight','bold')
                       .style('font-size',self.legend_text_size)



            
        },

        reDraw(){
            this.reDrawTree();
            this.reDrawScatter();
        }



    },


}
</script>

<style scoped lang="less">
    .corr-explorer-container{
        width: 100%;
        height: 580px;
        display: flex;
        flex-direction: row;
        align-items: center;
        .tree-container{
            flex:1 1 0;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            position:relative;
            margin:10px;

            .tree-title-container{
                border:2px solid white;
                border-radius: 10px;
                background-color :black;
                padding: 10px;
                margin:0 0 10px 0;
                b{
                    color:white;
                    font-size: 18px;
                }          
            }

            .tree-plot{
                width:700px;
                height:500px;
                border:2px solid lightgray;
                border-radius: 10px;
                background-color :white;

            }
        }

        .scatter-conatiner{
            flex:600px;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            margin:10px;
            .scatter-title-container{
                border:2px solid white;
                border-radius: 10px;
                background-color :black;
                padding: 10px;
                margin:0 0 10px 0;
                b{
                    color:white;
                    font-size: 18px;
                }          
            }

            .scatter-plot{
                height:500px;
                width:100%;
                border:2px solid lightgray;
                border-radius: 10px;
                background-color :white;
            }
        }
    }

    /deep/ .unselectNode{
        display: none;
    }

    /deep/ .el-dialog__body{//流水线对话框的主体（除去标题）
        background-color:#F5F5F5;
        border-top: 2px solid lightgray;
        padding: 10px 10px;
    }
    /deep/ .el-dialog__header span{
        font-size:20px;
    }

</style>