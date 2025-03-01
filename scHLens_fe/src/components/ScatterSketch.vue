<template>
    <div v-show="isShow" ref='sketch-container' class='sketch-container' :style="position">
        <svg ref="sketch-plot" class="sketch-plot"></svg>
    </div>

</template>

<script>
import Vue, { onBeforeMount } from "vue";
import * as d3 from "d3";
import {nanoid} from 'nanoid'
export default {
    name: "ScatterSketch",
    
    data(){
        return{
            nid:null,
            isShow:false,
            position:{
                left:'0px',
                top:'0px'
            },
            size:[0,0],
        }
    },

    methods:{
        
        draw(data,range=undefined,radius=undefined){//绘图函数
        // data:[{'x':...,'y':...,'color':...('#000000'),},...]
        // range:optional，data的range，对于图像的缩放性质有帮助，默认自适应：{'x':[xMin,xMax],'y':[yMin,yMax]}
        // radius:optional，圆点的尺寸大小，默认2
            //清空
            const svg = d3.select(this.$refs['sketch-plot'])
            
            svg.selectAll('*').remove();
            //判断特例
            if(data === undefined || data === null || data.length == 0)
                return;
            
            //配置图大小
            const width = this.size[0]
            const height = this.size[1]

            const padding={
                'top':10,
                'left':10,
                'bottom':10,
                'right':10,
            }

            //计算范围
            if(range === undefined){//未传入
                range = {
                    'x':[
                        Math.min(...data.map(v=>v.x)),
                        Math.max(...data.map(v=>v.x)),
                    ],
                    'y':[
                        Math.min(...data.map(v=>v.y)),
                        Math.max(...data.map(v=>v.y)),
                    ]
                }
            }
            //计算radius
            if(radius === undefined){//未传入
                // radius = 2;
                radius = Math.min(2000.0 / data.length,3)
            }

            //计算scale
            const xScale = d3.scaleLinear()
                .domain([...range.x])
                .range([radius + padding.left,width-radius-padding.right])
            const yScale = d3.scaleLinear()
                .domain([...range.y])
                .range([radius + padding.top,height-radius-padding.bottom])
            
            //绘图
            const plot = svg.append('g')
                .selectAll('*')
                .data(data)
                .join('circle')
                .attr('cx',d=>xScale(d.x))
                .attr('cy',d=>yScale(d.y))
                .attr('fill',d=>d.color)
                .attr('r',radius)
            

        },

        setSize(width,height){  
            d3.select(this.$refs['sketch-container']).style('width',`${width}px`).style('height',`${height}px`);
            this.size=[width,height]
        },

        setPos(left,top){
            this.position.left = `${left}px`
            this.position.top = `${top}px`
        },

        show(){
            this.isShow = true
        },

        hidden(){
            this.isShow = false
        }

    },

}
</script>

<style style scoped lang="less">
    .sketch-container{
        z-index: 99999;
        position: absolute;
        background-color: white;
        border: 2px solid rgb(200, 200, 200);
        box-shadow: 0px 0px 4px lightgray;
        border-radius: 2px;
        display: flex;
        justify-content: center;
        align-items: stretch;

        .sketch-plot{
            flex: 1 1; 
        }
    }
</style>