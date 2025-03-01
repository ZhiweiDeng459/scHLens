<template>
    <div class = "density-container">
        <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />
        <div class="resolution-slider">
                <span style="color:gray;font-size:16px;margin-right:5px;">Resolution:</span>
                <el-slider v-model="resolution" :min="2" :max="200" style="width:200px" @change="handleResolutionSliderChange"></el-slider>
        </div>
        <svg class="gene-density" ref="geneDensity" style="background-color: white"></svg>
    </div>
</template>

<script>
import * as d3 from "d3";
import SelfContextMenu from "@/components/SelfContextMenu"
import {requestGeneValueRange, requestGeneValueList} from '@/utils/interface'
import { Loading } from 'element-ui'
import {saveSvgAsPng} from 'save-svg-png-ext'
import eventBus from "@/utils/eventBus.js"
import { rgb } from "d3";

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
        infoPanel(){
            return this.$store.state.infoPanel
        },

    },
    data(){
        return {
            curGeneRange : [0,0],
            curGeneExpression : {},
            GeneColorScale : { //基因颜色转换器
                start_color:d3.rgb(211, 211, 211),
                mid_color:d3.rgb(255,255,0),
                end_color:d3.rgb(255,0,0),
                self:this,
                getColor:function(value){
                    if(value < this.self.curGeneRange[0]){
                        // return d3.rgb(255,255,255)
                        let tempScale = d3.scaleLinear().domain([this.self.curGeneRange[0]-1, 0]).range([0, 1]);
                        return d3.interpolate(d3.rgb(255,255,255),this.start_color)(tempScale(value))
                    }
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
            resolution:100,
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
            const svg = d3.select(self.$refs['geneDensity'])

            //清空
            svg.selectAll("*").remove();
            const legend = svg.append('g')
            const scatter = svg.append('g');

            /**计算初始数值 */

            const width = this.$refs['geneDensity'].clientWidth;
            const height = this.$refs['geneDensity'].clientHeight;

            const scatterWidth = width * 0.9;
            const scatterHeight = height;
            const scatterPadding = 60;

            const legendWidth = width * 0.1;
            const legendHeight = height;
            const legendPadding = 20;

            legend
                .attr("transform", `translate(${scatterWidth},${0})`)

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


            /**排除 */

            
            /**准备绘图数据 */
            //config
            //const density = Math.max(Math.min(150,Math.ceil(Math.sqrt(pointArr.length) * 3)),10) //10 - 100之间
            let density = self.resolution
            let widthResolution = density;//一行有多少个方块/有多少列 //TODO 自动设置resolution
            // const heightResolution = Math.ceil(1.0 * density * scatterHeight / scatterWidth) //TODO 自动设置resolution
            let heightResolution = density//一列有多少个方块/有多少行
            let unitWidth = 1.0 * scatterWidth / widthResolution
            let unitHeight = 1.0 * scatterHeight / heightResolution
            let contourData = new Array(widthResolution * heightResolution).fill(self.curGeneRange[0]-1);
            let contourData_count = new Array(widthResolution * heightResolution).fill(0);
            for(let point of pointArr){
                let i = Math.floor(posXScale(point.pos[0]) / unitWidth)
                let j = Math.floor(posYScale(point.pos[1]) / unitHeight)
                let index = j * widthResolution + i
                // // important !! 决定了单元格的值是最大值还是平均值
                // if(contourData[index] < self.curGeneExpression[point.id]){//最大值
                //     contourData[index] = self.curGeneExpression[point.id]
                // }
                if(contourData_count[index] == 0){
                    contourData[index] = self.curGeneExpression[point.id]
                }else{
                    contourData[index] += self.curGeneExpression[point.id]//平均值
                }
                contourData_count[index] += 1;
            }
            for(let i = 0; i < contourData.length;i++){//平均值
                if(contourData_count[i] != 0){
                    contourData[i] = 1.0 * contourData[i] / contourData_count[i]
                }
            }


            //blur
            // function check_out_of_range(i,j){
            //     if(i < 0)return false;
            //     if(i>=widthResolution)return false;
            //     if(j < 0)return false;
            //     if(j>=heightResolution)return false;
            //     return true;
            // }

            // let extend_iter_num =  1;
            // for(let iter = 0;iter < extend_iter_num;iter++){
            //     let newContourData = new Array(widthResolution * heightResolution).fill(self.curGeneRange[0]-1);
            //     let extend_size = 1;
            //     for(let i = 0; i < widthResolution;i++){
            //         for(let j = 0;j < heightResolution;j++){
            //             let index = j * widthResolution + i
            //             let value = contourData[index];
            //             let count = 1;
            //             for(let _i = i - extend_size; _i <= i+extend_size ;_i++){
            //                 for(let _j = j - extend_size; _j <= j+extend_size ;_j++){
            //                     if(_i == i && _j == j)continue;
            //                     if(!check_out_of_range(_i,_j))continue;
            //                     let dist = Math.sqrt(Math.pow(i - _i,2) + Math.pow(j - _j,2))
            //                     let weight = 1.0 * dist
            //                     value += contourData[_j * widthResolution + _i] * weight
            //                     count += weight
            //                 }
            //             }
            //             newContourData[index] = 1.0 * value / count;
            //         }
            //     }
            //     contourData = newContourData
            // }




            const path = d3.geoPath().projection(d3.geoIdentity().scale(unitWidth));
            const contours = d3.contours().size([widthResolution, heightResolution]).smooth(true);


            const color = function(value){
                return self.GeneColorScale.getColor(value)
            }

            let value_list = [self.curGeneRange[0]-1]
            for(let i = 0;i < 21;i++){
                value_list.push(self.curGeneRange[0] + i*0.05*(self.curGeneRange[1] - self.curGeneRange[0]))
            }


            svg.append("g")
                .selectAll()
                .data(value_list)
                .join("path")
                .attr("d", d => path(contours.contour(contourData, d)))
                .attr("fill", color)
                .classed('density-path',true)
                .on('mouseover',function(e,d){
                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'estimated value':d})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mousemove',function(e,d){
                    //show info
                    self.infoPanel.show()
                    self.infoPanel.setMessageData({'estimated value':d})
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                })
                .on('mouseout',function(){
                    //hidden info
                    self.infoPanel.hidden()

                })





            // /**绘图 */
            // let densityData = d3.contourDensity()
            //     .x(function(d) { return posXScale(d.pos[0]); })
            //     .y(function(d) { return posYScale(d.pos[1]); })
            //     .size([scatterWidth, scatterHeight])
            //     .bandwidth(20)(pointArr)

            // scatter
            //     .selectAll("path")
            //     .data(densityData)
            //     .enter().append("path")
            //     .attr("d", d3.geoPath())
            //     .attr("fill", function(d) { 
            //         return self.GeneColorScale.getColor(d.id); 
            //     })


            /**
             * 绘制legend
             */


            //清空
            legend.selectAll("*").remove();

            // //定义渐变色
            // const svgDefs = legendSVG.append("defs");
            // const testGradient = svgDefs.append("linearGradient")
            //     .attr("id", "densityGradient")
            //     .attr("x1",0)
            //     .attr("x2",0)
            //     .attr("y1",0)
            //     .attr("y2",1);
            // testGradient.append("stop").style("stop-color", "rgb(255, 0, 0)").attr("offset", "0");
            // testGradient.append("stop").style("stop-color", "rgb(255, 255, 0)").attr("offset", "0.5");
            // testGradient.append("stop").style("stop-color", "rgb(211, 211, 211)").style("opacity","0.4").attr("offset", "1.0");

            // //绘制表达值刻度
            // let ScaleX = 2,ScaleY = 10;
            // let ScaleWidth=10,ScaleHeigt = legendHeight - 2 * ScaleWidth;
            // let ScaleStrokeWidth = 1;
            // let ScaleLabelSize = 10;
            // legendSVG
            //     .append("rect")
            //     .style("fill","url(#densityGradient)")
            //     .attr("x", ScaleX)
            //     .attr("y", ScaleY)
            //     .attr("width", ScaleWidth)
            //     .attr("height", ScaleHeigt)
            //     .attr("stroke-width", ScaleStrokeWidth)
            //     .attr("stroke", "rgb(150, 150, 150)");
            
            // let NumScaleLine = 2;
            // let ScaleLineWidth = 5
            // for(let i = 0;i < NumScaleLine;i++){
            //     legendSVG.append("line")
            //         .attr("x1",ScaleX + ScaleWidth)
            //         .attr("y1",ScaleY + 0.5 * ScaleStrokeWidth + i * ((ScaleHeigt -  ScaleStrokeWidth)/(NumScaleLine-1)))
            //         .attr("x2",ScaleX + ScaleWidth + ScaleLineWidth)
            //         .attr("y2",ScaleY + 0.5 * ScaleStrokeWidth + i * ((ScaleHeigt -  ScaleStrokeWidth)/(NumScaleLine-1)))
            //         .style("stroke","rgb(150, 150, 150)")
            //         .style("stroke-width","1");
            //     legendSVG.append("text")
            //         .text(i==0?'max':'0')
            //         .attr("x",ScaleX + ScaleWidth + ScaleLineWidth + 1)
            //         .attr("y",ScaleY + i * (ScaleHeigt/(NumScaleLine-1)) + ScaleLabelSize * 0.4)
            //         .attr("font-size", ScaleLabelSize);
            // }
            //定义渐变色
            const svgDefs = legend.append("defs");
            const testGradient = svgDefs.append("linearGradient")
                .attr("id", "densityGradient")
                .attr("x1",0)
                .attr("x2",0)
                .attr("y1",0)
                .attr("y2",1);
            testGradient.append("stop").style("stop-color", self.GeneColorScale.end_color.formatHex()).attr("offset", "0");
            testGradient.append("stop").style("stop-color", self.GeneColorScale.mid_color.formatHex()).attr("offset", "0.5");
            testGradient.append("stop").style("stop-color", self.GeneColorScale.start_color.formatHex()).attr("offset", "1");


            //绘制表达值刻度
            let ScaleX = 2,ScaleY = 10;
            let ScaleWidth=10,ScaleHeigt = legendHeight - 2 * ScaleWidth;
            let ScaleStrokeWidth = 1;
            let ScaleLabelSize = 10;
            legend
                .append("rect")
                .attr("fill","url(#densityGradient)")
                .attr("x", ScaleX)
                .attr("y", ScaleY)
                .attr("width", ScaleWidth)
                .attr("height", ScaleHeigt)
                .attr("stroke-width", ScaleStrokeWidth)
                .attr("stroke", "rgb(150, 150, 150)");
            let NumScaleLine = 6;
            let ScaleLineWidth = 5
            
            //拿contourData的最大值作为legend的最大值

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
            //重绘所有视图
            eventBus.$emit("GeneProjectionViewRefreshingStart")
            if(this.curData.cellData.length == 0){//如果数据量为0
                d3.select(self.$refs['geneDensity']).selectAll('*').remove();
            }
            else{
                await this.updateCurGeneInfo();
                this.drawHeatMap();
            }

            eventBus.$emit("GeneProjectionViewRefreshingClose")

        },
        handleResolutionSliderChange(){
            this.drawHeatMap();
        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            
            //png保存
            // saveSvgAsPng(this.$refs.geneScatter, "geneScatter.png");

            //svg保存
            const svgDOM = this.$refs['geneDensity'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "geneDensity.svg";
            a.click();
            URL.revokeObjectURL(url)

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
.density-container{
    background-color: white;
    .resolution-slider{
        position:absolute;
        left:10px;
        top:0px;
        display: flex;
        align-items: center;
        width:300px;
    }

    .gene-density {
        width: 100%;
        height: 100%;
        /deep/ .density-path{
            // filter: blur(1px); /* 适度调小模糊值 */
        }

    }

}

</style>
