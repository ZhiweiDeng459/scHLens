<template>
    <div>
        <svg ref="TrajectoryInferenceCurve" id="TrajectoryInferenceCurve" style="background-color: white">
        </svg>
    </div>
</template>

<script>
import * as d3 from "d3";
import {saveSvgAsPng} from 'save-svg-png-ext'
export default {
    name:"TrajectoryInferenceCurve",
    props:['colorMode'],
    data(){
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
            this.drawScatter();
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
        draw(){
            //基础数据
            const width = this.$refs.TrajectoryInferenceCurve.clientWidth;
            const height = this.$refs.TrajectoryInferenceCurve.clientHeight;
            const padding = 50;
            let self = this;

            //初始化
            let svg = d3.select('#TrajectoryInferenceCurve')
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
            
            //绘制曲线
            const curve = svg.append('g')
            const line = d3.line().curve(d3.curveBasis).x(d=>posXScale(d[0])).y(d=>posYScale(d[1]))
            for(let l in this.TIData['shape']){
                let data = this.TIData['shape'][l]['curves'];
                curve.append('path')
                            .attr('fill','none')
                            .attr('stroke','black')
                            .attr('stroke-width',3)
                            .attr('d',line(data))
            }
            
            
        
        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            saveSvgAsPng(this.$refs.TrajectoryInferenceCurve, "TrajectoryInferenceCurve.png");            
        },
        reDraw(){
            //重绘所有
            this.draw();
        }
    },
    mounted(){
        this.reDraw();
    }
}
</script>

<style lang="less">
#TrajectoryInferenceCurve{
    height:100%;
    width:100%;
}
</style>