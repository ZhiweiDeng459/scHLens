<template>
    <div class = "density-container">
        <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />
        <div class="density">
            <div class="density-scatter" ref="heatmapContainer"></div>
            <svg class="density-scatter-legend" ref="heatmapLegend" style="background-color: white"></svg>
        </div>
    </div>
</template>

<script>
import * as d3 from "d3";
import h337 from "heatmap.js"
import SelfContextMenu from "@/components/SelfContextMenu"
import {requestGeneValueRange, requestGeneValueList} from '@/utils/interface'
import { Loading } from 'element-ui'
import {saveSvgAsPng} from 'save-svg-png-ext'
import eventBus from "@/utils/eventBus.js"

export default {
    name: "DensityScatter",
    components:{
        'SelfContextMenu':SelfContextMenu
    },
    computed: {
        curGeneName(){
            return this.$store.state.curGeneName
        },
        curData(){
            return this.$store.state.curData;
        },
        cellData(){
            return this.curData.cellData;
        },
        JobId(){
            return this.$store.state.JobId;
        },
        repaintTag(){
            return this.$store.state.repaintTag;
        },
    },
    data(){
        return {
            curGeneRange : [0,0],
            curGeneExpression : {},
            heatmapInstance:null,
            heatvalueConverter:{
                //避免value=0的时候，heatmap.js绘图机制出错
                //所有基因表达值范围强制变换到min-max
                min:1,
                max:100,
                self:this,
                convert:function(value){
                    return this.min + (value - this.self.curGeneRange[0])*(this.max - this.min)/(this.self.curGeneRange[1] - this.self.curGeneRange[0]);
                },
                reconvert:function(value){
                    return this.self.curGeneRange[0] + (value - this.min)*(this.self.curGeneRange[1] - this.self.curGeneRange[0])/(this.max - this.min);
                },
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

        'repaintTag.DensityScatter':{
            handler(){
                this.reDraw();
            }
        },
    },
    methods:{
        async updateCurGeneInfo(){
            //更新当前的基因信息：名称，范围，细胞表达值
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
        drawHeatMap(){
            const pointArr = JSON.parse(JSON.stringify(this.cellData));//深复制
            const self = this;
            /** 
            //排序保证绘制顺序
            pointArr.sort(function(cell1,cell2){
                if(self.curGeneExpression[cell1.id] < self.curGeneExpression[cell2.id])
                    return -1;
                else{
                    return 1;
                }
            });
            });**/ // TODO
            const width = this.$refs.heatmapContainer.clientWidth;
            const height = this.$refs.heatmapContainer.clientHeight;

            const scatterWidth = width * 0.9;
            const scatterHeight = height;
            const scatterPadding = 60;

            // let minX = Math.min(...pointArr.map(d=>d.pos[0]));
            // let maxX = Math.max(...pointArr.map(d=>d.pos[0]));
            // let minY = Math.min(...pointArr.map(d=>d.pos[1]));
            // let maxY = Math.max(...pointArr.map(d=>d.pos[1]));
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
            
            //组装热力图的点
            let points = [];
            for(let i = 0; i < pointArr.length;i++){
                let point = {
                    x:posXScale(pointArr[i].pos[0]),
                    y:posYScale(pointArr[i].pos[1]),
                    value:this.heatvalueConverter.convert(this.curGeneExpression[pointArr[i].id]),
                };
                points.push(point);
            }
            let data = {
                max:this.heatvalueConverter.convert(this.curGeneRange[1]),
                min:this.heatvalueConverter.convert(this.curGeneRange[0]),
                data:points,
            };
            this.heatmapInstance.setData(data);
        },
        drawHeatMapLegend(){
            const height = this.$refs.heatmapLegend.clientHeight;
            const legendSVG = d3.select(`.density-scatter-legend`);

            //清空
            legendSVG.selectAll("*").remove();

            //定义渐变色
            const svgDefs = legendSVG.append("defs");
            const testGradient = svgDefs.append("linearGradient")
                .attr("id", "densityGradient")
                .attr("x1",0)
                .attr("x2",0)
                .attr("y1",0)
                .attr("y2",1);
            testGradient.append("stop").style("stop-color", "rgb(255, 0, 0)").attr("offset", "0");
            testGradient.append("stop").style("stop-color", "rgb(255, 255, 0)").attr("offset", "0.5");
            testGradient.append("stop").style("stop-color", "rgb(211, 211, 211)").style("opacity","0.4").attr("offset", "1.0");

            //绘制表达值刻度
            let ScaleX = 2,ScaleY = 10;
            let ScaleWidth=10,ScaleHeigt = height - 2 * ScaleWidth;
            let ScaleStrokeWidth = 1;
            let ScaleLabelSize = 10;
            legendSVG
                .append("rect")
                .style("fill","url(#densityGradient)")
                .attr("x", ScaleX)
                .attr("y", ScaleY)
                .attr("width", ScaleWidth)
                .attr("height", ScaleHeigt)
                .attr("stroke-width", ScaleStrokeWidth)
                .attr("stroke", "rgb(150, 150, 150)");
            
            let NumScaleLine = 2;
            let ScaleLineWidth = 5
            for(let i = 0;i < NumScaleLine;i++){
                legendSVG.append("line")
                    .attr("x1",ScaleX + ScaleWidth)
                    .attr("y1",ScaleY + 0.5 * ScaleStrokeWidth + i * ((ScaleHeigt -  ScaleStrokeWidth)/(NumScaleLine-1)))
                    .attr("x2",ScaleX + ScaleWidth + ScaleLineWidth)
                    .attr("y2",ScaleY + 0.5 * ScaleStrokeWidth + i * ((ScaleHeigt -  ScaleStrokeWidth)/(NumScaleLine-1)))
                    .style("stroke","rgb(150, 150, 150)")
                    .style("stroke-width","1");
                legendSVG.append("text")
                    .text(i==0?'max':'0')
                    .attr("x",ScaleX + ScaleWidth + ScaleLineWidth + 1)
                    .attr("y",ScaleY + i * (ScaleHeigt/(NumScaleLine-1)) + ScaleLabelSize * 0.4)
                    .attr("font-size", ScaleLabelSize);
            }
        },
        async reDraw(){
            //重绘所有视图
            eventBus.$emit("GeneProjectionViewRefreshingStart")
            if(this.curData.cellData.length == 0){//如果数据量为0
                d3.select('.density-scatter-legend').selectAll('*').remove();
                d3.select('.density-scatter').selectAll('*').remove();
            }
            else{
                await this.updateCurGeneInfo();
                if(this.heatmapInstance === null){
                    let gradient = {
                        0.0 :"rgb(211, 211, 211)",
                        0.5 :"rgb(255,255,0)",
                        1.0 :"rgb(255,0,0)",                 
                    };
                    this.heatmapInstance = h337.create({
                        container: this.$refs.heatmapContainer,
                        radius:5,
                        gradient:gradient,
                        minOpacity: 0.4,
                        backgroundColor:'white',
                    });
                }
                this.drawHeatMap();
                this.drawHeatMapLegend();
            }

            eventBus.$emit("GeneProjectionViewRefreshingClose")

        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            // saveSvgAsPng(this.$refs.heatmapContainer, "GeneHeatMap.png");
            const canvasElement = d3.select(".heatmap-canvas").node();
            const MIME_TYPE = "image/png";
            const imgURL = canvasElement.toDataURL(MIME_TYPE);
            const dlLink = document.createElement('a');
            dlLink.setAttribute("download","GeneHeatMap.png")
            dlLink.href = imgURL;
            dlLink.dataset.downloadurl = [MIME_TYPE,dlLink.download,dlLink.href].join(':')
            document.body.appendChild(dlLink)
            dlLink.click();
            document.body.removeChild(dlLink)

        },
        menuMounted(){

        },
    },
    mounted(){
        this.reDraw();
    }
};
</script>

<style scoped lang="less">
.density{
    width:100%;
    height:100%;
    background-color:white;
    .density-scatter {
        float:left;
        width: 88%;
        height: 100%;
    }

    .density-scatter-legend {
        float:right;
        width:10%;
        height: 100%;
    }
}
</style>
