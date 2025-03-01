<template>
  <div>
    <svg ref="CellChatEdgePlot" id="CellChatEdgePlot" style="background-color: white">
    </svg>
  </div>
</template>

<script>
import * as d3 from "d3";
import {saveSvgAsPng} from 'save-svg-png-ext'
export default {
    name:'CellChatEdgePlot',

    data() {
        return {
            mode:'count',

        }
    },
    computed:{
        curData() {
            return this.$store.state.curData;
        },
        groups(){
            return this.curData.groups;
        },
        CC(){
            return this.curData.CC;
        }
    },
    watch:{
        curData(){
            if (this.curData === undefined || this.curData === null) return;
            this.reDraw();
        },

    },
    methods:{
        drawPlot(){
            const width = this.$refs.CellChatEdgePlot.clientWidth;
            const height = this.$refs.CellChatEdgePlot.clientHeight;
            const padding = 15;
            
            const size = Math.min(width - 2 *padding,height - 2 * padding)
            let self = this;

            const centerX = width * 0.5;
            const centerY = height * 0.5;

            
            const svg = d3.select('#CellChatEdgePlot')
            svg.selectAll('*').remove();


            //计算边界
            let cellCount = 0; //细胞总数
            this.groups.forEach(element => {
                cellCount += element['size']
            });
            
            const maxGroupSize = Math.max(...this.groups.map(v=>v['size'])) //最大的group的细胞数
 
            const maxDistance = 0.5 * size / (1 + Math.sin(Math.PI / this.groups.length)); //保证聚类圆不相交，不溢出的又可接受的聚类圆心到点的距离，TODO 这里设置为默认距离

            const maxRadius = 0.5 *  maxDistance * Math.sin(Math.PI / this.groups.length); //在maxDistance的前提下，可接受的最大聚类半径。如果设定距离小于maxDistance，那么半径也要适当缩短

            //整理数据
            let rawData = this.CC[this.mode];
            let names =  this.groups.map(item=>item['id'])
            let attachGroupData = this.groups.map((item,i)=>{ //groups的辅助信息
                return {
                    'id':item['id'],
                    'angle': i * 2 * Math.PI / self.groups.length,//圆心所在的角度
                    'x' : centerX + Math.cos( i * 2 * Math.PI / self.groups.length ) * maxDistance,//圆心的横坐标
                    'y' : centerY + Math.sin( i * 2 * Math.PI / self.groups.length ) * maxDistance, //圆心的纵坐标
                    'r' : 1.0 * this.groups[i].size / maxGroupSize * maxRadius, //圆的半径
                }
            })
            let drawData = []
            let maxWeight = 0;
            for(let source of names){
                for(let target of names){
                    maxWeight = Math.max(maxWeight,rawData[source][target])
                    drawData.push({
                        'source':source,
                        'target':target,
                        'value':rawData[source][target]
                    })
                }
            }


            //绘制连接线
            const link = svg.append('g')
            link.selectAll('path')
                .data(drawData)
                .join('path')
                .attr('d',function(d,i){
                    //计算端点的位置
                    let source = attachGroupData.find(v=>v.id == d.source)
                    let target = attachGroupData.find(v=>v.id == d.target)
                    let distance = Math.sqrt( 
                        Math.pow(source.x - target.x,2) +  Math.pow(source.y - target.y,2)
                        )

                    let _source = { //位于两个圆心连接线上，在source圆边上的点
                        'x':source.x + (1.0 * source.r / distance) * (target.x - source.x),
                        'y':source.y + (1.0 * source.r / distance) * (target.y - source.y),
                    }
                    let _target = { //位于两个圆心连接线上，在source圆边上的点
                        'x':target.x + (1.0 * target.r / distance) * (source.x - target.x),
                        'y':target.y + (1.0 * target.r / distance) * (source.y - target.y),
                    }
                    
                    let _k = 0.2 //_center向圆心O偏移的度数
                    let _center = { //连接线的中段控制点，用于产生平滑的弯曲
                        'x':_k * centerX + (1 - _k) / 2.0 * (source.x + target.x),
                        'y':_k * centerY + (1 - _k) / 2.0 * (source.y + target.y)
                    }
                    const link = d3.line()
                                .curve(d3.curveBasis)
                                .x(d=>d.x)
                                .y(d=>d.y)

                    return link([source,_center,target])

                })
                .attr('stroke',function(d){
                    return self.groups.find(v=>v.id == d.source).color
                })
                .attr('stroke-width','10px')
                .attr('opacity',0.5)
                .attr('fill',function(d){
                    return 'none'
                })


            // //整理ribbon生成器
            // const ribbon = d3.ribbonArrow()
            //                 .source(d=>{
            //                     return {
            //                         'angle':attachGroupData.find(v=>v.id == d.source).angle,
            //                         'ribbonR':maxDistance - attachGroupData.find(v=>v.id == d.source).r,
            //                         'value':d.value
            //                     }
            //                 })
            //                 .target(d=>{
            //                     return {
            //                         'angle':attachGroupData.find(v=>v.id == d.target).angle,
            //                         'ribbonR':maxDistance - attachGroupData.find(v=>v.id == d.target).r,
            //                         'value':d.value
            //                     }
            //                 })
            //                 .radius(d=>{
            //                     return d.ribbonR
            //                 })
            //                 .startAngle(d=>{
            //                     return d.angle + 0.5 * Math.PI - d.value / maxWeight * 0.08
            //                 })
            //                 .endAngle(d=>{
            //                     return d.angle + 0.5 * Math.PI + d.value / maxWeight * 0.08
            //                 })      

            // //绘制ribbon
            // const arrow = svg.append('g')
            //                 .attr("transform", `translate(${0.5*width},${0.5*height})`)
            // arrow.selectAll('path')
            //     .data(drawData)
            //     .join('path')
            //     .attr('d',function(d){
            //         //设定ribbon属性
            //         return ribbon(d);             
            //     })
            //     // .attr('stroke','white')
            //     // .attr('stroke-width','0.8px')
            //     .attr('opacity',0.5)
            //     .attr('fill',d=>this.groups.find(v=>v.id==d.source).color)

        
            //绘制圆
            const node = svg.append('g');
            node.selectAll('circle')
                .data(this.groups)
                .join('circle')
                .attr('cx',function(d,i){
                    return attachGroupData[i].x;
                })
                .attr('cy',function(d,i){
                    return attachGroupData[i].y;
                })
                .attr('r',function(d,i){
                    return attachGroupData[i].r
                })
                .attr('fill',d=>d.color)


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
#CellChatEdgePlot{
    height:100%;
    width:100%;

}

</style>