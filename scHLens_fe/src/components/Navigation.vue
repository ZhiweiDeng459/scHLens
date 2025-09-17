<template>
    <div ref="TreeContainer" class="TreeContainer">
        <svg ref="HierarchyPlot" class="HierarchyPlot"></svg>
        <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
            />
        <ScatterSketch ref="scatter-sketch"/>
    </div>

</template>

<script>
import Vue from "vue";
import { Button ,Image,Notification,MessageBox} from "element-ui";
import * as d3 from "d3";
import curNodeImg from '@/assets/icons/curNode.svg'
import NodeImg from '@/assets/icons/Node.svg'
import RedDeleteImg from '@/assets/icons/RedDelete.svg'
import BlackDeleteImg from '@/assets/icons/BlackDelete.svg'

import SelfContextMenu from "@/components/SelfContextMenu"
import ScatterSketch from "@/components/ScatterSketch"
import eventBus from "@/utils/eventBus.js"


import {requestDeleteView} from "@/utils/interface"

Vue.component(Button.name, Button);
Vue.component(Image.name, Image);
Vue.component(Notification.name,Notification)
Vue.component(MessageBox.name,MessageBox);



export default {
    name: "Navigation",
    components:{
        SelfContextMenu,
        ScatterSketch,
    },
    props:[
        'callMergeOptions'
    ],
    data() {
        return {
            menuItems:[
                {
                    'name':'Merge',
                    'icon':'icons/merge.svg',
                    'callback':this.callMergeOptions
                },
                {
                    'name':'Select Cells in this Node',
                    'icon':'icons/selectWithNode.svg',
                    'callback':this.selectWithNode
                }
            ],
            sketchWidth:100,
            sketchHeight:100,
            nodeSize:20,
            SelectRingRadius:10,
            setMergeRoot:undefined,//重设mergeRoot的函数
        };
    },
    computed: {
        curData(){
            return this.$store.state.curData;
        },
        ViewIdList() {
            return this.$store.state.dataList.map((item) => item.ViewId);
        },
        dataList(){
            return this.$store.state.dataList;
        },
        projectTree() {
            return this.$store.state.projectTree;
        },
        chosenViews(){
            return this.$store.state.chosenViews;
        },
        mergeRoot(){
            return this.$store.state.mergeRoot;
        },
        repaintTag(){
            return this.$store.state.repaintTag;
        },
        JobId(){
            return this.$store.state.JobId;
        }

    },
    watch: {
        projectTree(){
            this.reDraw()
        },
        curData(){
            this.reDraw()
        },
        'repaintTag.ProjectionTree':{
            handler(){
                this.reDraw();
            }
        },
    },

    methods: {
        //绘制
        drawTree() {
            if(this.projectTree['root'] === null)
                return ;
            const svg = d3.select(".HierarchyPlot"); 


            if(svg === undefined)
                return;
            if(this.$refs.HierarchyPlot === undefined)
                return;


            svg.selectAll('*').remove();


            const padding = {
                'top' : 40,
                'left' : 30,
                'right' : 30,
                'bottom' : 50,
            }

            const height =  this.$refs.HierarchyPlot.clientHeight;
            const width = this.$refs.HierarchyPlot.clientWidth;


            const tree = d3.tree()
                .size([width - padding.left - padding.right,height - padding.top - padding.bottom])


            let self = this;

            //准备数据
            const rootNode = this.projectTree['root'];
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
                root.y = 0.5 * (height - padding.top - padding.bottom)
            }


            self.nodeSize = Math.min(30,0.3 * (height - padding.top - padding.left) / (root.height + 1))
            self.SelectRingRadius = 0.5 * self.nodeSize

            //select Ring
            const selectRing = svg.append('g')
            root.each(d => {
                selectRing
                .append('circle')
                .attr('r',self.SelectRingRadius)
                .attr('cx',d.x + padding.left)
                .attr('cy',d.y + padding.top)
                .attr('stroke',function(){
                    if(self.mergeRoot == d.data.name){
                        return '#409EFF'
                    }
                    else{
                        return '#32CD32'
                    }
                })
                .attr('stroke-width',4)
                .attr('fill',()=>{
                    return 'white'
                })
                .classed('unselectNode',function(){
                    let findResult = self.chosenViews.find((item)=>{
                            return item == d['data']['name']
                        })
                    if(findResult === undefined)
                        return true;
                    else{
                        return false;
                    }
                })
            });


            function setMergeRoot(){//根据当前的chosenViews设置mergeRoot
                ////找到最高节点
                let minDepth = 9999;
                let tempMergeRoot = undefined;
                let tempMergeRootNode = undefined;
                for(let v of self.chosenViews){
                    let item = root.find(node=>{
                        return node.data.name == v
                    }) 
                    if(item.depth < minDepth){
                        minDepth = item.depth
                        tempMergeRoot = item.data.name
                        tempMergeRootNode = item
                    }
                }
                ////校验是否具有合并的可能性
                if(tempMergeRoot !== undefined){
                    for(let v of self.chosenViews){
                        if(tempMergeRootNode.find(node=>node.data.name == v) === undefined){
                            tempMergeRoot = undefined;
                            break;
                        }
                    }
                    if(self.chosenViews.length < 2){
                        tempMergeRoot = undefined;
                    }
                    self.$store.state.mergeRoot = tempMergeRoot
                    self.$store.commit('updateMergeRoot',tempMergeRoot)
                }

            }
            self.setMergeRoot = setMergeRoot;

            //entity
            const Entity = svg.append('g')
            root.each(d=>{
                Entity.append('image')
                    .attr('width',self.nodeSize)
                    .attr('height',self.nodeSize)
                    .attr('x',function(){
                        return d.x + padding.left - this.getBoundingClientRect().width * 0.5
                    })
                    .attr('y',function(){
                        return d.y + padding.top - this.getBoundingClientRect().height * 0.5;
                    })
                    .attr('href',function(){
                        return (d['data']['name'] == self.curData['ViewId'] ?  curNodeImg : NodeImg)
                    })
                    .style('cursor','pointer')
                    .on('click',function(){//选择/取消选择元素
                        let findIndex = self.chosenViews.indexOf(d['data']['name'])
                        if(findIndex !== -1){ //删除选中的视图

                            self.chosenViews.splice(findIndex,1);
                        }
                        else{ //选中该视图
                            self.chosenViews.push(d['data']['name']);
                        }

                        //寻找mergeRoot
                        self.setMergeRoot();

                        self.reDraw();
                    })
                    .on("contextmenu", function (event) {//切换视图
                        event.preventDefault();
                        let ViewId = d['data']['name']
                        self.$store.commit("toggleCurData", ViewId);
                        //对于无数据的节点，要展示提示
                        if(self.dataList.find(v=>v.ViewId == ViewId).cellData.length == 0){
                            Notification({
                            title: 'Notification',
                            dangerouslyUseHTMLString: true,
                            message: `There are no cells existing in this node <b>${self.curData.ViewId}</b>, so its views won't show anything.`,
                            type: 'warning',
                            duration: 15000
                            });
                        }
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

                        const containerWidth = self.$refs['TreeContainer'].getBoundingClientRect().width
                        const containerHeight = self.$refs['TreeContainer'].getBoundingClientRect().height
                        let pos = [event.offsetX + 10,event.offsetY + 10]
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
                        const containerWidth = self.$refs['TreeContainer'].getBoundingClientRect().width
                        const containerHeight = self.$refs['TreeContainer'].getBoundingClientRect().height
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
                        .style('font-size','12px')
                        .attr('x',function(){
                            return d.x + padding.left - this.getBoundingClientRect().width * 0.5
                        })
                        .attr('y',function(){
                            return d.y + padding.top + self.SelectRingRadius + this.getBoundingClientRect().height
                        })
                        .style('font-family','YaHei')
                        .style('font-weight','bold')
            })


            //arrow
            let arrowSize = 13
            svg.append('defs')
               .append('marker')
               .attr('id','navigation_arrow')
               .attr('markerWidth',arrowSize)
               .attr('markerHeight',arrowSize)
               .attr('refX',0)
               .attr('refY',0.5 * arrowSize)
               .attr("markerUnits","userSpaceOnUse")
               .attr('orient', 'auto')
               .append('path')
               .attr('fill', '#434343')
               .attr('d', `M 0,0 L ${arrowSize},${0.5*arrowSize} L 0,${arrowSize}`)


            //link
            const links = svg.append('g')
                .selectAll('path')
                .data(root.links())
                .join('path')
                .attr('fill','none')
                .attr('stroke','#666666')
                .attr('stroke-width',5)
                .attr('marker-end','url(#navigation_arrow)')
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
                                            'y':parseFloat(sourceText.getAttribute('y')) + 5
                                        }
                                   })
                                   .target(function(d){
                                        
                                        return {
                                            'x':d.target.x,
                                            'y':d.target.y - self.SelectRingRadius + padding.top - arrowSize - 1,
                                        }
                                   })
                                   .x(d=>d.x + padding.left)
                                   .y(d=>d.y)
                    return link(data)
                })


            //delete button
            const DeleteButton = svg.append('g')
            root.each(d=>{
                DeleteButton.append('image')
                    .attr('width',10)
                    .attr('height',10)
                    .attr('x',function(){
                        return d.x + padding.left + 9
                    })
                    .attr('y',function(){
                        return d.y + padding.top - this.getBoundingClientRect().height -9;
                    })
                    .attr('href',function(){
                        return BlackDeleteImg
                    })
                    .style('cursor','pointer')
                    .on("mouseover",function(event){
                        d3.select(this).attr('href',RedDeleteImg)
                    })
                    .on("mouseout",function(event){
                        d3.select(this).attr('href',BlackDeleteImg)
                    })
                    .on("click",function(event){//删除该节点，及其下属节点
                        self.deleteView(d['data']['name'])
                    })

            })

        },

        deleteView(ViewId){//删除视图及其下属视图
            const self = this;
            //如果是根节点，那么禁止删除，并且提示
            if(ViewId == self.projectTree['root']['id']){
                this.$message({
                    'message':'The root node is not allowed to be deleted',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            self.$store.commit('deleteDataObj',ViewId);
            //重新设置mergeRoot
            self.setMergeRoot();
            //刷新视图
            self.reDraw();
            //服务器同步删除
            requestDeleteView(self.JobId,ViewId)
                .then((res)=>{
                    console.log(res)
                })
                .catch((err)=>{
                    console.log(err)
                    this.$message({
                        'message':'The operation to delete the view failed to upload to the server',
                        'type':'error',
                        'showClose':true,
                    })
                })


        },

        menuMounted(_this,root,parent){

        },
        reDraw(){
            // eventBus.$emit("ProjectionTreeRefreshingStart")
            this.drawTree();
            // this.$refs['TreeContainer'].hidden()
            // eventBus.$emit("ProjectionTreeRefreshingClose")
        },
        toggleView(ViewId) {
            this.$store.commit("toggleCurData", ViewId);
        },
        selectWithNode(){
            //如果选择了多个节点或者没有选择节点的话，报错
            if(this.chosenViews.length == 0){
                // this.$message({
                //     'message':'Please select at least one node to perform this operation',
                //     'type':'error',
                //     'showClose':true,
                // })
                MessageBox.alert(
                            `<strong>Details：</strong>Please select at least one node to perform this operation`,
                            'Error: Select Cells in this Node',
                            {
                                dangerouslyUseHTMLString: true,
                                type: 'error',
                            }
                        );
                return;
            }
            //获取选择节点
            for(let targetViewId of this.chosenViews){
                let targetData = this.dataList.find(v=>{
                    return v['ViewId'] == targetViewId
                });
                //获取CellIds
                let CellIds = targetData.cellData.map(v=>v.id)
                //发送选择请求
                eventBus.$emit('selectCellsById',CellIds);
            }
        }

    },

    mounted(){
    }
};
</script>


<style scoped lang="less">
.TreeContainer{
    position: relative;
    .HierarchyPlot{
        height: 100%;
        width: 100%;
    }
}



/deep/ .unselectNode{
    display: none;
}



</style>
