<template>
  <div>
    <svg ref="TrajectoryInferencePlot" id="TrajectoryInferencePlot" style="background-color: white">
    </svg>
  </div>
</template>

<script>
import * as d3 from "d3";
import {saveSvgAsPng} from 'save-svg-png-ext'
export default {
    name:'TrajectoryInferencePlot',
    data() {
        return {

        }
    },
    computed:{
        curData() {
            return this.$store.state.curData;
        },
        TIData(){
            return this.curData.TI;
        },
        groups(){
            return this.curData.groups;
        }
    },
    watch:{
        curData(){
            if (this.curData === undefined || this.curData === null) return;
            this.drawPlot();
        },
        groups:{
            //主要是监控组名被修改
            deep: true,
            handler(newValue,oldValue) {
                if (oldValue === "undefined" || oldValue === "null") return;
                this.reDraw();
            },
        },
    },
    methods:{
        drawPlot(){
            const width = this.$refs.TrajectoryInferencePlot.clientWidth;
            const height = this.$refs.TrajectoryInferencePlot.clientHeight;

            const padding = 50;
            
            const maxRadius = 20;
            const minRadius = 5;
            
            const maxEdgeWidth = 10;

            let self = this;

            let svg = d3.select('#TrajectoryInferencePlot')
            svg.selectAll("*").remove();


            const scatterData = JSON.parse(JSON.stringify(this.TIData.scatter));
            const meanArr = JSON.parse(JSON.stringify(this.TIData.mean))
            const edges = []
            let TD = JSON.parse(JSON.stringify(this.TIData.connectivities));
            for(let i = 0;i < meanArr.length;i++){
                for(let j = i + 1; j < meanArr.length; j++){
                    if(TD[i][j] != 0)
                        edges.push({
                            'node1': i,//TODO 这里的节点不是id，而是与meanArr顺序绑定
                            'node2': j,
                            'weight': TD[i][j]
                        })
                }
            }

            let maxX = Math.max(...scatterData.map(item=>{
                return item['pos'][0]
            }));
            let minX = Math.min(...scatterData.map(item=>{
                return item['pos'][0]
            }));
            let maxY = Math.max(...scatterData.map(item=>{
                return item['pos'][1]
            }));
            let minY = Math.min(...scatterData.map(item=>{
                return item['pos'][1]
            }));
            let maxSize = Math.max(...this.groups.map(item=>{
                return item['size']
            }))
            let minSize = Math.min(...this.groups.map(item=>{
                return item['size']
            }))

            const posXScale = d3
                .scaleLinear()
                .domain([minX, maxX])
                .range([padding, width - padding]);
            const posYScale = d3
                .scaleLinear()
                .domain([minY, maxY])
                .range([padding, height - padding]);
            
            const RadiusScale = d3
                .scaleLinear()
                .domain([minSize,maxSize])
                .range([minRadius,maxRadius])

            //绘制links
            const links = svg.append('g')
            links.selectAll("line")
                .data(edges)
                .enter()
                .append("line")
                .attr("x1",(d) => posXScale(meanArr[d['node1']]['X']))
                .attr("y1",(d) => posYScale(meanArr[d['node1']]['Y']))
                .attr("x2",(d) => posXScale(meanArr[d['node2']]['X']))
                .attr("y2",(d) => posYScale(meanArr[d['node2']]['Y']))
                .classed("ti-plot-line",true)
                .style("stroke-width",(d) => maxEdgeWidth * d['weight'])       
            
            //绘制nodes
            const nodes = svg.append('g')
            nodes.selectAll("circle")
                .data(meanArr)
                .enter()
                .append("circle")
                .attr("cx", (d) => {
                    return posXScale(d['X'])
                })
                .attr("cy", (d) => posYScale(d['Y']))
                .attr("r", (d,i)=>{
                    return RadiusScale(self.curData.groups[i].size)
                })
                .attr("fill", (d,i) => this.groups[i].color);
            //绘制text
            const texts = svg.append('g')
            texts.selectAll("text")
                .data(JSON.parse(JSON.stringify(this.groups)))
                .enter()
                .append("text")
                .text(d=>d.name)
                .classed("ti-plot-text",true)
                .attr("x",function(d,i){
                    return posXScale(meanArr[i]['X']) - this.getBoundingClientRect().width * 0.52;
                })
                .attr("y",function(d,i){
                    return posYScale(meanArr[i]['Y']) + this.getBoundingClientRect().height * 0.35;
                })

        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            saveSvgAsPng(this.$refs.TrajectoryInferencePlot, "TrajectoryInferencePlot.png");            
        },
        reDraw(){
            //重绘所有
            this.drawPlot();
        }
    },
    mounted(){
        this.reDraw()
    }
}
</script>

<style lang="less">
#TrajectoryInferencePlot{
    height:100%;
    width:100%;
    .ti-plot-line{
        stroke: gray;
    }
    .ti-plot-text{
        font-size: 20px;
        font-weight: 700;
    }   
}

</style>