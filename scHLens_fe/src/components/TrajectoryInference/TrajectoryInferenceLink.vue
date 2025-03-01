<template>
  <div>
    <svg ref="TrajectoryInferenceLink" id="TrajectoryInferenceLink" style="background-color: white">
    </svg>
  </div>
</template>

<script>
import * as d3 from "d3";
import {saveSvgAsPng} from 'save-svg-png-ext'
export default {
    name:'TrajectoryInferenceLink',
    props:['colorMode'],
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
        scatterData(){
            return this.curData.cellData
        },
        colorScheme(){
            return this.curData.colorScheme
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
        colorMode(){
            this.reDraw();
        }
    },
    methods:{
        drawPlot(){
                      //基础数据
            const width = this.$refs.TrajectoryInferenceLink.clientWidth;
            const height = this.$refs.TrajectoryInferenceLink.clientHeight;
            const padding = 50;
            let self = this;

            //初始化
            let svg = d3.select('#TrajectoryInferenceLink')
            svg.selectAll("*").remove();
            
            //绘图数据
            const scatterData = JSON.parse(JSON.stringify(this.scatterData));


            //比例尺
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
            const posXScale = d3
                .scaleLinear()
                .domain([minX, maxX])

                .range([padding, width - padding]);
            const posYScale = d3
                .scaleLinear()
                .domain([minY, maxY])
                .range([padding, height - padding]);
            
            //绘制点
            const scatter = svg.append('g')
            scatter.selectAll("circle")
                .data(scatterData)
                .enter()
                .append("circle")
                .attr("cx", (d) => posXScale(d['pos'][0]))
                .attr("cy", (d) => posYScale(d['pos'][1]))
                .attr("r", 3)
                .attr("fill", function(d,i){
                    if(self.colorMode == 'group')
                        return self.groups.find(e=>d.group == e.id).color
                    else if(self.colorMode == 'pseudotime')
                        return self.TIData.PseudotimeColor[i]
                });
    
            //绘制中心连线
            const link = svg.append('g')
            link.selectAll('path')
            const line = d3.line()
                .x(d=>{
                    let group = self.groups.find(v=>v.id==d)
                    return posXScale(group['centerX'])
                })
                .y(d=>{
                    let group = self.groups.find(v=>v.id==d)
                    return posYScale(group['centerY'])
                })
            for(let l in this.TIData['shape']){
                let data = this.TIData['shape'][l]['lineages'];
                link.append('path')
                            .attr('fill','none')
                            .attr('stroke','black')
                            .attr('stroke-width',3)
                            .attr('d',line(data))
            }
            //绘制中心点
            const innerRadius = 10
            const outerRadius = 15;
            const outerCenter = svg.append('g')
            outerCenter.selectAll('circle')
                .data(this.groups)
                .join('circle')
                .attr("cx",(d) => posXScale(d['centerX']))
                .attr("cy",(d) => posYScale(d['centerY']))
                .attr("r",outerRadius)
                .attr("fill", 'white')
                .attr('stroke','black')
                .attr('stroke-width',2 + 'px')
            const innerCenter = svg.append('g')
            innerCenter.selectAll('circle')
                .data(this.groups)
                .join('circle')
                .attr("cx",(d) => posXScale(d['centerX']))
                .attr("cy",(d) => posYScale(d['centerY']))
                .attr("r",innerRadius)
                .attr("fill", (d) => d.color)
                .attr('stroke','black')
                .attr('stroke-width',2 + 'px')
                .attr('fill-opacity',0.6)

        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            saveSvgAsPng(this.$refs.TrajectoryInferenceLink, "TrajectoryInferenceLink.png");            
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
#TrajectoryInferenceLink{
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