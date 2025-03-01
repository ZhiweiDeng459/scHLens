<template>
  <div>
    <svg ref="ChordPlot" id="ChordPlot" style="background-color: white">
    </svg>
  </div>
</template>

<script>
import * as d3 from "d3";
import {chord} from "d3-chord"
import {saveSvgAsPng} from 'save-svg-png-ext'
export default {
    name:'CellChatChordPlot',

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
            const width = this.$refs.ChordPlot.clientWidth,
            height = this.$refs.ChordPlot.clientHeight,
            padding = 50;
            let self = this;
            
            
            //整理数据
            let rawData = this.CC[this.mode];
            
            let names =  this.groups.map(item=>item['id'])

            let matrix = []
            for(let i = 0;i < names.length;i++){
                let temp = []
                for(let j = 0;j < names.length;j++){
                    temp.push(rawData[names[i]][names[j]])
                }
                matrix.push(temp)
            }



            let innerRadius = Math.min(width, height) * 0.5 - 90;
            let outerRadius = innerRadius + 10
            let ArrowRootRaduis = 10
            let arcOuter = d3.arc()
                .innerRadius(innerRadius + ArrowRootRaduis + 2)
                .outerRadius(outerRadius + ArrowRootRaduis + 2)
            let ribbon = d3.ribbonArrow()
                .radius(innerRadius - 1)
                .padAngle(1 / innerRadius)
            let chord = d3.chordDirected()
                .padAngle(10 / innerRadius)
                .sortSubgroups(d3.descending)
                .sortChords(d3.descending)


            //开始绘图
            let svg = d3.select('#ChordPlot')
            svg.selectAll("*").remove();
            
            const chords = chord(matrix);


            const rootSVG = svg.append('g')
                .attr("transform", `translate(${0.5*width},${0.5*height})`)
                  

            const group = rootSVG.append("g")
                .attr("font-size", 10)
                .attr("font-family", "sans-serif")
                .selectAll("g")
                .data(chords.groups)
                .join("g");
            
            const arrowRoot = rootSVG.append("g")
                .selectAll("g")
                .data(chords)
                .join("g")

            group.append("path")
                .attr("fill", d => {return this.groups[d.index]['color']})
                .attr("d", arcOuter);

            arrowRoot.append("path")
                .attr("fill",d=>this.groups[d.target.index]['color'])
                .attr("d",function(d,i){
                    const rootArc = d3.arc()
                        .innerRadius(innerRadius)
                        .outerRadius(outerRadius)
                        .startAngle(d.source.startAngle)
                        .endAngle(d.source.endAngle)
                    return rootArc();
                })


            group.append("text")
                .each(d => (d.angle = (d.startAngle + d.endAngle) / 2))
                .attr("dy", "0.35em")
                .attr("transform", d => `
                    rotate(${(d.angle * 180 / Math.PI - 90)})
                    translate(${outerRadius + 15})
                    ${d.angle > Math.PI ? "rotate(180)" : ""}
                `)
                .attr("text-anchor", d => d.angle > Math.PI ? "end" : null)
                .text(d => names[d.index])
                .style("font-size","22px")

            group.append("title")
                .text(d => `${names[d.index]}
            ${d3.sum(chords, c => (c.source.index === d.index) * c.source.value)} outgoing →
            ${d3.sum(chords, c => (c.target.index === d.index) * c.source.value)} incoming ←`);

            rootSVG.append("g")
                .attr("fill-opacity", 0.75)
                .selectAll("path")
                .data(chords)
                .join("path")
                .style("mix-blend-mode", "multiply")
                .attr("fill", d => this.groups[d.target.index]['color'])
                .attr("d", ribbon)
                .append("title")
                .text(d => `${names[d.source.index]} → ${names[d.target.index]} ${d.source.value}`);
                
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
#ChordPlot{
    height:100%;
    width:100%;

}

</style>