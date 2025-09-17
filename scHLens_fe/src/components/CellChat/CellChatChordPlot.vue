<template>
  <div ref="cellchat-chord-container" class="cellchat-chord-container">
    <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />
    <svg ref="ChordPlot" id="ChordPlot" style="background-color: white">
    </svg>
  </div>
</template>

<script>
import * as d3 from "d3";
import {chord} from "d3-chord"
import {saveSvgAsPng} from 'save-svg-png-ext'
import eventBus from "@/utils/eventBus.js"
import SelfContextMenu from "@/components/SelfContextMenu"

export default {
    name:'CellChatChordPlot',
    props:['setInteractionTable'],
    components:{
        'SelfContextMenu':SelfContextMenu
    },

    data() {
        return {
            mode:'weight',
            menuItems:[
                {
                    'name':'Save this Image',
                    'icon':'icons/save_as_image.svg',
                    'callback':()=>{
                        this.saveToFile();
                    }
                }
            ],


        }
    },
    computed:{
        curData() {
            return this.$store.state.curData;
        },
        groups(){
            return this.curData.groups;
        },
        curGeneName(){
            return this.$store.state.curGeneName
        },
        CC(){
            return this.curData.CC;
        },
        repaintTag(){
            return this.$store.state.repaintTag;
        },
        infoPanel(){
            return this.$store.state.infoPanel
        },
    },
    watch:{
        curData(){
            if (this.curData === undefined || this.curData === null) return;
            this.reDraw();
        },
        'repaintTag.CellChatEdge':{
            handler(){
                //TODO更新分数
                this.reDraw();
            }
        },
        groups: {
            //监视
            deep: true,
            handler() {
                if (this.groups === "undefined" || this.groups === "null") return;
                this.reDraw()
            },
        },
    },
    methods:{
        drawPlot(){

            //计算大小
            const raw_height = this.$refs['cellchat-chord-container'].clientHeight
            const raw_width = this.$refs['cellchat-chord-container'].clientWidth


            const width = raw_width;
            const height = raw_height;

            let svg = d3.select('#ChordPlot')
            //设置视图大小
            svg.attr("width",width)
            svg.attr('height',height)

            const padding = 20;
            const self = this;
            
            
            //整理数据
            let rawData = this.CC[this.mode];
            
            let names =  this.groups.map(item=>item['id'])
            let annos = this.groups.map(item=>item['name'])

            let matrix = []
            for(let i = 0;i < names.length;i++){
                let temp = []
                for(let j = 0;j < names.length;j++){
                    temp.push(rawData[names[i]][names[j]])
                }
                matrix.push(temp)
            }


            let interaction_map = {}
            for(let i = 0;i < names.length;i++){
                interaction_map[names[i]] = {} 
                for(let j = 0;j < names.length;j++){
                    interaction_map[names[i]][names[j]] = []
                }
            }
            for(let i = 0;i < self.CC['source'].length;i++){
                let source = self.CC['source'][i]
                let target = self.CC['target'][i]
                let ligand = self.CC['ligand'][i]
                let receptor = self.CC['receptor'][i]
                let score = self.CC['prob'][i]
                interaction_map[source][target].push({
                    'ligand':ligand,
                    'receptor':receptor,
                    'score':score.toPrecision(4)
                })
            }


            let innerRadius = Math.min(width, height) * 0.5 - 60;
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
            svg.selectAll("*").remove();
            


            const chords = chord(matrix);

            const zoomSVG = svg.append('g')
                               .attr("class", "cellchat-chord-zoomLayer")
            
            //定义zoom行为
            const zoom = d3.zoom()
                               .scaleExtent([0.2, 20])  // 设置缩放范围，最小 0.5 倍，最大 5 倍
                               .on("zoom", (event) => {
                                    zoomSVG.attr("transform", event.transform);  // 平移缩放
                                });
            svg.call(zoom)

            const rootSVG = zoomSVG.append('g')
                .attr("transform", `translate(${0.5*width},${0.5*height})`)
                  
            //外圈
            const group = rootSVG.append("g")
                .attr("font-size", 10)
                .attr("font-family", "sans-serif")
                .selectAll("g")
                .data(chords.groups)
                .join("g");

            group.append("path")
                .attr("fill", d => {return this.groups[d.index]['color']})
                .attr("d", arcOuter);

            //内圈
            const arrowRoot = rootSVG.append("g")
                .selectAll("g")
                .data(chords)
                .join("g")

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
                .text(d => annos[d.index])
                .style("font-size","14px")
                .style('font-family','YaHei')
                .style('font-weight','bold')

            //arrow
            rootSVG.append("g")
                .attr("fill-opacity", 0.75)
                .selectAll("path")
                .data(chords)
                .join("path")
                .style("mix-blend-mode", "multiply")
                .attr("fill", d => this.groups[d.target.index]['color'])
                .attr("d", ribbon)
                .classed('cell-chat-chord-arrow',true)
                .on("mouseover",function(e,d){
                    //突出显示
                    d3.select(this)
                        .attr("stroke","black")
                        .attr("stroke-width",2)
                    //显示信息
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({
                        'Source':self.groups[d.source.index].name,
                        'Target':self.groups[d.target.index].name,
                        'Strength':matrix[d.source.index][d.target.index].toPrecision(3),
                    })
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 30)
                })
                .on("mousemove",function(e,d){
                    d3.select(this)
                        .attr("stroke","black")
                        .attr("stroke-width",2)
                    //显示信息
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({
                        'Source':self.groups[d.source.index].name,
                        'Target':self.groups[d.target.index].name,
                        'Strength':matrix[d.source.index][d.target.index].toPrecision(3),
                    })
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 30)

                })
                .on("mouseout",function(e,d){
                    d3.select(this)
                        .attr("stroke",null)
                        .attr("stroke-width",null)
                    self.infoPanel.hidden()

                })
                .on("click",function(e,d){
                    let table_data = interaction_map[self.groups[d.source.index].id][self.groups[d.target.index].id]

                    self.setInteractionTable(table_data)
                })


        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            // //png保存
            // saveSvgAsPng(this.$refs.dotplot, "dotplot.png");
            //svg保存
            const svgDOM = this.$refs['ChordPlot'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "Cell Chat View - ChordPlot.svg";
            a.click();
            URL.revokeObjectURL(url)
        },
        reDraw(){
            eventBus.$emit('CellChatViewRefreshingStart');
            //重绘所有
            this.drawPlot();
            eventBus.$emit('CellChatViewRefreshingClose');
        },
        menuMounted(){

        },
    },
    mounted(){
        this.reDraw()
    }
}
</script>

<style lang="less">
#ChordPlot{
    display: block;
    cursor:move;
}
.cell-chat-chord-arrow{
    cursor: pointer;
}
.cellchat-chord-container{
    display: inline-block;
    display: flex;
    position: relative;
    align-items: stretch;

}


</style>