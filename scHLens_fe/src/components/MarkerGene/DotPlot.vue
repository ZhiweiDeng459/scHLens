<template>
    <div>
        <div class="dotplot-container">
            <el-scrollbar ref="dotplot-scroll" class="dotplot-main-view">
                <svg class="dotplot" ref="dotplot" style="background-color: white"></svg>
            </el-scrollbar>
           <SelfContextMenu
                :items = menuItems
           />
        </div>
    </div>
</template>

<script>
import * as d3 from "d3";
import Vue from "vue";
import { Scrollbar,Loading } from 'element-ui';
import {saveSvgAsPng} from 'save-svg-png-ext'
import SelfContextMenu from "@/components/SelfContextMenu"
import eventBus from "@/utils/eventBus.js"

Vue.use(Scrollbar);
export default {
    name: "DotPlot",
    components:{
        SelfContextMenu
    },
    data() {
        return {
            rankRange: [1, 10],
            menuItems:[
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
            return this.$store.state.curData
        },
        cellData() {
            return this.curData.cellData;
        },
        groups() {
            return this.curData.groups;
        },
        marker(){
            return this.curData.MK;
        },
        dendrogram(){
            return this.curData.dendrogram;
        },
        repaintTag(){
            return this.$store.state.repaintTag;
        },
        infoPanel(){
            return this.$store.state.infoPanel
        },
        curGeneName(){
            return this.$store.state.curGeneName;
        },

    },
    watch: {
        marker() {
            if (this.marker === "undefined" || this.marker === "null") return;
            this.reDraw();
        },
        groups: {
            //监视group的名称在scatter部分被修改
            deep: true,
            handler() {
                if (this.groups === "undefined" || this.groups === "null") return;
                this.reDraw();
            },
        },
        'repaintTag.DotPlot':{
            handler(){
                this.reDraw();
            }
        },
    },
    methods: {


        drawDotPlotHorizontal(){//水平方向作图
            //绘图先决条件判断
            if(this.marker === undefined || this.marker === null || this.marker.length == 0)
                return ;



            const self = this;
            const data = this.marker;  //TODO 返回的marker基因存在重复的情况
            
            const marker_names = data.map(v=>v.name)
                    
            const group_names = this.groups.map(v=>v.name)
            const group_ids = this.groups.map(v=>v.id)
        
            const svg = d3.select(`.dotplot`);
            svg.selectAll("*").remove();

            const axisMargin = { top: 0, right: 40, bottom: 30, left: 15 }
            const cellSize = 35                                                                                           
            const radius = 14
                
            const xAxisLength = (marker_names.length - 1) * cellSize;
            const yAxisLength = (this.groups.length - 1) * cellSize;


            let maxMean = -Infinity;
            let minMean = Infinity;
            for(let v of data.map((item) => item.means).flat()){
                if(v >= maxMean)
                    maxMean = v;
                if(v <= minMean)
                    minMean = v;
            }
            
            const colorScale = d3.scaleLinear().domain([minMean, maxMean]).range([0, 0.8]);

            // const colorScale = function(value){
            //     let valueScale = d3.scaleLinear().domain([minMean,maxMean]).range([0, 1]);
            //     let index = valueScale(value)
            //     if(value > 0.5){
            //         return d3.interpolate(d3.rgb(69, 180, 194),d3.rgb(255, 248, 37))((index-0.5)*2);  
            //     }
            //     else{
            //         return d3.interpolate(d3.rgb(8, 29, 88),d3.rgb(69, 180, 194))(index*2);
            //     }
            // }

            //draw gene axis
            const geneScale = d3.scalePoint()
                                .domain(marker_names)
                                .range([0,xAxisLength])
            const geneAxis = d3.axisTop(geneScale)
                                .ticks(marker_names.length)
                                .tickFormat(function(d){
                                    let gene_name = d.substring(0,d.lastIndexOf("@"))
                                    let tick_name = gene_name
                                    if(tick_name.length > 9){
                                        tick_name = tick_name.substring(0,9) + '...'
                                    }
                                    return tick_name
                                })
            svg.append('g')
                .classed('geneAxis',true)
                .call(geneAxis)
            
            //draw group axis
            const groupScale = d3.scalePoint()
                                .domain(group_ids)
                                .range([0,yAxisLength])
            const groupAxis = d3.axisLeft(groupScale).ticks(group_ids.length).tickFormat((d)=>{
                                    return this.groups.find(g=>g.id==d).name
                                })
            svg.append('g')
                .classed('groupAxis',true)
                .call(groupAxis)

            //draw the group domain
            const groupInterval = svg.append('g')
            const groupIntervalLine = groupInterval.append('g')
                                                    .selectAll('g')
                                                    .data(this.groups)
                                                    .join('g')
            const groupIntervalLabel = groupInterval.append('g')
                                                    .selectAll('g')
                                                    .data(this.groups)
                                                    .join('g')
            groupIntervalLine
                .append('line')
                .attr("x1",(d,i)=>geneScale(marker_names[i * marker_names.length / this.groups.length]))
                .attr("y1", 35)
                .attr("x2",(d,i)=>geneScale(marker_names[(i + 1) * marker_names.length / this.groups.length - 1]))
                .attr("y2", 35)
                .attr("stroke-width", 2)
                .attr("stroke", "#434343")
            groupIntervalLine
                .append('line')
                .attr("x1",(d,i)=>geneScale(marker_names[i * marker_names.length / this.groups.length]))
                .attr("y1",35)
                .attr("x2",(d,i)=>geneScale(marker_names[i * marker_names.length / this.groups.length]))
                .attr("y2",45)
                .attr("stroke-width", 2)
                .attr("stroke", "#434343")
            groupIntervalLine
                .append('line')
                .attr("x1",(d,i)=>geneScale(marker_names[(i + 1) * marker_names.length / this.groups.length - 1]))
                .attr("y1",35)
                .attr("x2",(d,i)=>geneScale(marker_names[(i + 1) * marker_names.length / this.groups.length - 1]))
                .attr("y2",45)
                .attr("stroke-width", 2)
                .attr("stroke", "#434343")
            groupIntervalLabel //group标签
                .append('text')
                .text(d=>d.name)
                .attr("x",(d,i)=>0.5 * (geneScale(marker_names[i * marker_names.length / this.groups.length]) + geneScale(marker_names[(i + 1) * marker_names.length / this.groups.length - 1])))
                .attr("y",25)
                .attr('text-anchor',"middle")
                .attr('font-size','13px')
                .style("font-family","YaHei")
                .style("font-weight","bold")



            //set tick size and rotate
            d3.select('.geneAxis') //gene axis
                .selectAll('g')
                .selectAll('text')
                .attr('font-size','12px')
                .style("transform","rotate(35deg)")
                .style("transform",function(){
                    return `translate(${0},${-this.getBoundingClientRect().height * 0.4}px)` + d3.select(this).style('transform')
                })
                .style("font-family","YaHei")
                .style("font-weight","bold")
                .on("mouseover",function(e,d){ //在悬浮时突出显示
                    d3.select(this)
                        .style("fill","#3b80ee")
                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'gene name':d.substring(0,d.lastIndexOf("@"))})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)
                })
                .on("mousemove",function(e,d){ //在悬浮时突出显示
                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'gene name':d.substring(0,d.lastIndexOf("@"))})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)
                })
                .on("mouseout",function(){
                    d3.select(this)
                        .style("fill","black")
                    //hidden info
                    self.infoPanel.hidden()
                        
                })
            d3.select('.groupAxis') //group axis
                .selectAll('g')
                .selectAll('text')
                .attr('font-size','14px')
                .style("font-family","YaHei")
                .style("font-weight","bold")


            //move element
            const xMove = svg.select('.groupAxis').node().getBoundingClientRect().width + axisMargin.left;
            const yMove = svg.select('.geneAxis').node().getBoundingClientRect().height + groupInterval.node().getBoundingClientRect().height + axisMargin.top;
            svg.select('.geneAxis').attr('transform',`translate(${xMove + 0.5 * cellSize},${yMove})`)
            svg.select('.groupAxis').attr('transform',`translate(${xMove},${yMove + 0.5 * cellSize})`)
            groupInterval.attr('transform',`translate(${xMove + 0.5 * cellSize},0)`)

            //draw rect
            svg.append("rect")
                .attr("x",xMove)
                .attr("y",yMove)
                .attr("width",xAxisLength + cellSize)
                .attr("height",yAxisLength + cellSize)
                .attr("stroke-width", 2)
                .attr("stroke", "#434343")
                .attr("fill", "rgba(255, 255, 255, 0)");

            //draw dot and dynamic legend
            const dot = svg.append('g')
                        .attr('transform',`translate(${0.5 * cellSize},${0.5 * cellSize})`)
            
            const dyanmic_legend = svg.append('g')
                            .attr('transform',`translate(${0.5 * cellSize},${0.5 * cellSize})`) 
            
            //dot cell
            const dot_cell = dot.selectAll('*')
                .data(marker_names)
                .join('g')
                .selectAll('*')
                .data((marker,i)=>data[i].means.map((d,j)=>{
                    return {
                        'gene':marker,
                        'mean':d,
                        'fraction':data[i].fraction[j],
                        'group_index':j
                    }
                }))
                .join('circle')
                .attr("cx",(d,i)=>geneScale(d.gene) + xMove)
                .attr("cy",(d,i)=>groupScale(this.groups[i].id) + yMove)
                .attr("fill",(d,i)=>d3.interpolateYlOrRd(colorScale(d.mean)))
                .attr("r",(d,i)=>{
                    return radius * d.fraction
                })


            //dynamic legend
            const legend_g = dyanmic_legend.selectAll('*')
                .data(marker_names)
                .join('g')
                .selectAll('*')
                .data((marker,i)=>data[i].means.map((d,j)=>{
                    return {
                        'gene':marker,
                        'mean':d,
                        'fraction':data[i].fraction[j],
                        'group_index':j
                    }
                }))
                .join('g')
                .attr('opacity',0)
                .on('mouseover',function(e,d){
                    d3.select(this).attr('opacity',1)

                    let gene_name = d['gene'].substring(0,d['gene'].lastIndexOf("@"))
                    let group_name = self.groups[d.group_index]['name']

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'group':group_name,'marker gene':gene_name,'mean':d['mean'].toFixed(3),'fraction':d['fraction'].toFixed(3)})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mousemove',function(e,d){
                    d3.select(this).attr('opacity',1)

                    let gene_name = d['gene'].substring(0,d['gene'].lastIndexOf("@"))
                    let group_name = self.groups[d.group_index]['name']

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'group':group_name,'marker gene':gene_name,'mean':d['mean'].toFixed(3),'fraction':d['fraction'].toFixed(3)})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mouseout',function(){
                    d3.select(this).attr('opacity',0)
                    //hidden info
                    self.infoPanel.hidden()

                })


            legend_g.append('text')
                .text((d)=>(d.fraction * 100).toFixed(0) + '%')
                .style('font-family','YaHei')
                .style('font-weight','bold')
                .style('font-size','9px')
                .attr('text-anchor','middle')
                .attr('dominant-baseline','central')
                .attr("x",function(d){
                    return geneScale(d.gene) + xMove //- this.getBoundingClientRect().width * 0.5
                })
                .attr("y",function(d){
                    return groupScale(self.groups[d.group_index].id) + yMove// + this.getBoundingClientRect().height * 0.5
                })


            legend_g.append('circle')
                .attr("cx",(d,i)=>geneScale(d.gene) + xMove)
                .attr("cy",(d,i)=>groupScale(this.groups[i].id) + yMove)
                .attr("fill-opacity",0)
                .attr("r",(d,i)=>{
                    return radius * 1
                })
                .style('stroke','black')
                .style('stroke-width',2)


            //set size
            const width = xMove + xAxisLength + cellSize + axisMargin.right; 
            const height = yMove + yAxisLength + cellSize + axisMargin.bottom;
            svg.attr("width", width);
            svg.attr("height", height);


            //set event
            d3.select('.geneAxis')
                .selectAll('g')
                .style("cursor", "pointer")
                .on("click", function (event,d) {
                    //裁剪出基因真正名称
                    let real_gene = d.substring(0,d.lastIndexOf("@"))

                    if(self.curGeneName.includes(real_gene)){//基因已经有了该基因
                        self.$message({
                            'message':'This Gene has been added',
                            'type':'error',
                            'showClose':true,
                        })
                        return;
                    }

                    //上传基因修改事件
                    self.$store.commit("addToCurGeneName", real_gene);
                    
                });

        },


        reDraw(){

            eventBus.$emit("MarkerGeneViewRefreshingStart")
            if(this.curData.cellData.length == 0 || this.marker.length == 0){//如果数据量为0
                d3.select('.dotplot').selectAll('*').remove();
            }
            else{
                this.drawDotPlotHorizontal()
            }

            this.$refs['dotplot-scroll'].update();
            
            eventBus.$emit("MarkerGeneViewRefreshingClose")

        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            // //png保存
            // saveSvgAsPng(this.$refs.dotplot, "dotplot.png");
            //svg保存
            const svgDOM = this.$refs['dotplot'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "Marker Gene View.svg";
            a.click();
            URL.revokeObjectURL(url)
        },
    },
    mounted(){
        this.reDraw();
    }
};
</script>

<style lang="less">
.dotplot-container{
    display: flex;
    position: relative;
    align-items: stretch;
    height: 100%;
    width: 100%;
    flex-direction: column;
    .dotplot-main-view {
        flex: 1 1;
        //overflow-x: auto;
        //overflow-y: auto;
        border-right: 2px solid rgb(200, 200, 200);
        min-height: 0;
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


    .dotplot-legend {
        flex:0 200px;
        //border: 2px solid rgb(200, 200, 200);
        .stop-left {
            stop-color: rgb(255, 245, 240);
        }
        .stop-middle {
            stop-color: rgb(249, 105, 76);
        }
        .stop-right {
            stop-color: rgb(103, 0, 13);
        }
        .gradient-red {
            fill: url(#mainGradient);
        }
    }
}
</style>
