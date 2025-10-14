<template>
    <div class = "gene-scatter-container">
        <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />
        <svg class="gene-scatter" ref="geneScatter" style="background-color: white"></svg>
    </div>
</template>

<script>
import * as d3 from "d3";
import SelfContextMenu from "@/components/SelfContextMenu"
import {requestGeneValueRange, requestGeneValueList} from '@/utils/interface'
import {saveSvgAsPng} from 'save-svg-png-ext'
import eventBus from "@/utils/eventBus.js"

export default {
    name: "GeneScatter",

    computed: {
        curGeneName(){
            return this.$store.state.curGeneName;
        },
        curData(){
            return this.$store.state.curData;
        },
        dataset(){
            return this.curData.paramsObj['dataset'];
        },
        cellData(){
            return this.curData.cellData;
        },
        JobId(){
            return this.$store.state.JobId
        },
        repaintTag(){
            return this.$store.state.repaintTag;
        },
        infoPanel(){
            return this.$store.state.infoPanel
        },
        groups(){
            return this.curData.groups;
        },

    },
    data() {
        return {
            curGeneRange : [0,0],
            curGeneExpression : {},
            GeneColorScale : {
                start_color:d3.rgb(8, 29, 88),
                mid_color:d3.rgb(69, 180, 194),
                end_color:d3.rgb(255, 248, 37),
                self:this,
                getColor:function(value){
                    //全部表达值相同的的情况
                    if(this.self.curGeneRange[1] == this.self.curGeneRange[0]){
                        return this.start_color;
                    }
                    //表达值不同的情况
                    let GeneScale = d3.scaleLinear().domain([this.self.curGeneRange[0], this.self.curGeneRange[1]]).range([0, 1]);
                    let index = GeneScale(value);
                    if(index > 0.5){
                        return d3.interpolate(this.mid_color,this.end_color)((index-0.5)*2);
                    }
                    else{
                        return d3.interpolate(this.start_color,this.mid_color)(index*2);
                    }

                },
                getIndex:function(value){
                    let GeneScale = d3.scaleLinear().domain([this.self.curGeneRange[0], this.self.curGeneRange[1]]).range([0, 1]);
                    return GeneScale(value);
                },
                InvertToValue:function(index){
                    let GeneScale = d3.scaleLinear().domain([this.self.curGeneRange[0], this.self.curGeneRange[1]]).range([0, 1]);
                    return GeneScale.invert(index);
                }
            },
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
    components:{
        'SelfContextMenu':SelfContextMenu
    },
    watch: {
        async cellData() {
            if (this.cellData === undefined || this.cellData === null) return;
            this.reDraw();
        },
        curGeneName:{
            deep:true,
            async handler(){
                this.reDraw();
            }
        },

        'repaintTag.GeneScatter':{
            handler(){
                this.reDraw();
            }
        },
    },
    methods: {
        async updateCurGeneInfo(){
            //更新当前的基因信息：范围，细胞表达值
            //获取当前基因范围信息
            await requestGeneValueRange(this.JobId,this.curData.ViewId,this.curGeneName)
            .then((response) => {
                    this.curGeneRange = response.data;
                })
                .catch((err) => {
                    console.log(err);
                });
            //获取基因表达信息
            await requestGeneValueList(this.JobId,this.curData.ViewId,this.curGeneName)
            .then((response)=>{
                    this.curGeneExpression = response.data;
            })
            .catch((err) => {
                console.log(err);
            });
        },
        drawScatterCircle() {
            const pointArr = JSON.parse(JSON.stringify(this.cellData));//注意，这里是一次深copy
            const self = this;
            const svg = d3.select(`.gene-scatter`);

            //清空
            svg.selectAll("*").remove();
            const legend = svg.append('g');
            const scatter = svg.append('g');
            
            const height = this.$refs.geneScatter.clientHeight;   //svg总高度
            const width = height //svg总宽度

            const scatterWidth = width * 0.9;
            const scatterHeight = height;
            const scatterPadding = 60;

            const legendWidth = width * 0.1;
            const legendHeight = height;
            const legendPadding = 20;

            legend
                .attr("transform", `translate(${scatterWidth + 60},${0})`)
            
            /**
             * 
             * 绘制scatter
             * 
             */
            
            let minX = this.curData.raw_embedding_range['x'][0]
            let maxX = this.curData.raw_embedding_range['x'][1]
            let minY = this.curData.raw_embedding_range['y'][0]
            let maxY = this.curData.raw_embedding_range['y'][1]

            const posXScale = d3
                .scaleLinear()
                .domain([minX, maxX])
                .range([scatterPadding, scatterWidth -  scatterPadding]);
            const posYScale = d3
                .scaleLinear()
                .domain([minY, maxY])
                .range([scatterPadding, scatterHeight - scatterPadding]);
            
            //绘制基因表达图的点
            //按照基因表达程度从小到大排序，完成堆叠次序
            pointArr.sort(function(cell1,cell2){
                if(self.curGeneExpression[cell1.id] < self.curGeneExpression[cell2.id])
                    return -1;
                else{
                    return 1;
                }
            });
            scatter.selectAll("circle")
                .data(pointArr)
                .enter()
                .append("circle")
                .attr("cx", (d) => posXScale(d.pos[0]))
                .attr("cy", (d) => posYScale(d.pos[1]))
                .attr("fill", (d) => this.GeneColorScale.getColor(this.curGeneExpression[d.id]))
                .attr("r", Math.max(Math.min(3000.0 / pointArr.length,3),1))
                .attr("id", (d) => d.id)
                .attr("class", (d) => `label${d.group}`)
                .on('mouseover',function(e,d){

                    let anno = self.groups.find(g=>{
                        return g.id == d.group
                    }).name

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'id':d.id,'group':anno,'value':self.curGeneExpression[d.id]})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mousemove',function(e,d){
                    let anno = self.groups.find(g=>{
                        return g.id == d.group
                    }).name

                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'id':d.id,'group':anno,'value':self.curGeneExpression[d.id]})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mouseout',function(){
                    //hidden info
                    self.infoPanel.hidden()

                })


            /**
             * 
             * 绘制legend
             * 
             */
            //定义渐变色
            const svgDefs = legend.append("defs");
            const testGradient = svgDefs.append("linearGradient")
                .attr("id", "geneGradient")
                .attr("x1",0)
                .attr("x2",0)
                .attr("y1",0)
                .attr("y2",1);
            testGradient.append("stop").style("stop-color", self.GeneColorScale.end_color.formatHex()).attr("offset", "0");
            testGradient.append("stop").style("stop-color", self.GeneColorScale.mid_color.formatHex()).attr("offset", "0.5");
            testGradient.append("stop").style("stop-color", self.GeneColorScale.start_color.formatHex()).attr("offset", "1");


            //绘制表达值刻度
            let ScaleX = 2,ScaleY = 10;
            let ScaleWidth=10,ScaleHeigt = height - 2 * ScaleWidth;
            let ScaleStrokeWidth = 1;
            let ScaleLabelSize = 17;
            legend
                .append("rect")
                .attr("fill","url(#geneGradient)")
                .attr("x", ScaleX)
                .attr("y", ScaleY)
                .attr("width", ScaleWidth)
                .attr("height", ScaleHeigt)
                .attr("stroke-width", ScaleStrokeWidth)
                .attr("stroke", "rgb(150, 150, 150)");
            let NumScaleLine = 6;
            let ScaleLineWidth = 5
            for(let i = 0;i < NumScaleLine;i++){
                legend.append("line")
                    .attr("x1",ScaleX + ScaleWidth)
                    .attr("y1",ScaleY + 0.5 * ScaleStrokeWidth + i * ((ScaleHeigt -  ScaleStrokeWidth)/(NumScaleLine-1)))
                    .attr("x2",ScaleX + ScaleWidth + ScaleLineWidth)
                    .attr("y2",ScaleY + 0.5 * ScaleStrokeWidth + i * ((ScaleHeigt -  ScaleStrokeWidth)/(NumScaleLine-1)))
                    .style("stroke","rgb(150, 150, 150)")
                    .style("stroke-width","1");
                legend.append("text")
                    .text(((this.curGeneRange[1] - this.curGeneRange[0])*(1 - i*0.2) + this.curGeneRange[0]).toPrecision(2))
                    .attr("x",ScaleX + ScaleWidth + ScaleLineWidth + 1)
                    .attr("y",ScaleY + i * (ScaleHeigt/(NumScaleLine-1)) + ScaleLabelSize * 0.4)
                    .attr("font-size", ScaleLabelSize);
            }    

            
        },
        async reDraw(){
            /**
             * 全部重绘
             */
            eventBus.$emit("GeneProjectionViewRefreshingStart")
            
            if(this.curData.cellData.length == 0){//如果数据量为0
                d3.select('.gene-scatter').selectAll('*').remove();
            }
            else{
                await this.updateCurGeneInfo();
                // this.drawScatterLegend();
                this.drawScatterCircle();
            }
            eventBus.$emit("GeneProjectionViewRefreshingClose")
        },
        menuMounted(){

        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            
            //png保存
            // saveSvgAsPng(this.$refs.geneScatter, "geneScatter.png");

            //svg保存
            const svgDOM = this.$refs['geneScatter'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "Gene Projection View - Scatter.svg";
            a.click();
            URL.revokeObjectURL(url)

        }
    },
    mounted(){
        this.reDraw()
    }
};
</script>

<style lang="less">
.gene-scatter-container{
    background-color: white;
    .gene-scatter {
        width: 100%;
        height: 100%;
        .gene-left {
            stop-color: rgb(255, 248, 37);
        }
        .gene-middle {
            stop-color: rgb(69, 180, 194);
        }
        .gene-right {
            stop-color: rgb(8, 29, 88);
        }
    }
}
</style>
