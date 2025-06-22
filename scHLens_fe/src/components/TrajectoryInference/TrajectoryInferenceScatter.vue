<template>
    <div>
        <svg ref="TrajectoryInferenceScatter" id="TrajectoryInferenceScatter" style="background-color: white">
        </svg>
    </div>
</template>

<script>
import * as d3 from "d3";
import {saveSvgAsPng} from 'save-svg-png-ext'
export default {
    name:"TrajectoryInferenceScatter",
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
    },
    methods:{
        drawScatter(){
            const width = this.$refs.TrajectoryInferenceScatter.clientWidth,
            height = this.$refs.TrajectoryInferenceScatter.clientHeight,
            padding = 50;

            let self = this;

            let svg = d3.select('#TrajectoryInferenceScatter')
            svg.selectAll("*").remove();
            
            const scatterData = JSON.parse(JSON.stringify(this.TIData['scatter']));
            const meanArr = JSON.parse(JSON.stringify(this.TIData.mean))

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
            svg.selectAll("circle")
                .data(scatterData)
                .enter()
                .append("circle")
                .attr("cx", (d) => posXScale(d['pos'][0]))
                .attr("cy", (d) => posYScale(d['pos'][1]))
                .attr("r", 3)
                .attr("fill", (d) => this.groups.find(e=>d.group == e.id).color);
    
            //绘制text
            svg.selectAll("text")
                .data(JSON.parse(JSON.stringify(this.groups)))
                .enter()
                .append("text")
                .text(d=>d.name)
                .classed("ti-scatter-text",true)
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
            saveSvgAsPng(this.$refs.TrajectoryInferenceScatter, "TrajectoryInferenceScatter.png");            
        },
        reDraw(){
            //重绘所有
            this.drawScatter();
        }
    },
    mounted(){
        this.reDraw();
    }
}
</script>

<style lang="less">
#TrajectoryInferenceScatter{
    height:100%;
    width:100%;
    .ti-scatter-text{
        font-size: 20px;
        font-weight: 700;
    }   
}
</style>