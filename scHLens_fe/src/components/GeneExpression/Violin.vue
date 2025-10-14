<template>
    <div class="violin-container" ref="violinContainer">
        <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />
        <el-scrollbar ref="violin-scroll" class = "violin-plot-container" >
            <svg class="violin" ref="violin" style="background-color: white"></svg>
        </el-scrollbar>
    </div>
</template>

<script>
import {requestGeneValueRange, requestGeneValueList} from '@/utils/interface';
import { Scrollbar,Loading } from 'element-ui'
import {saveSvgAsPng} from 'save-svg-png-ext'
import * as d3 from "d3";
import Vue from "vue"
import SelfContextMenu from "@/components/SelfContextMenu"
import eventBus from "@/utils/eventBus.js"


Vue.use(Scrollbar);
export default {
    name: "Violin",
    props:['mode'],
    components:{
        SelfContextMenu
    },
    data(){
        return {
            //基本框
            width:0,
            height:0,
            unitWidth:0,
            unitPlotMaxWidth:0,
            margin:{},
            //数据相关
            curGeneRange : [0,0],
            curGeneExpression : {},
            curViolinData:[],
            curDensity:[],
            //绘图相关
            thresholds:[],
            wiggle:null,
            xScale:null,
            yScale:null,
            xAxis:null,
            yAxis:null,
            areal:null,
            arear:null,
            //菜单相关
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
    computed:{
        curData(){
            return this.$store.state.curData;
        },
        groups(){
            return this.curData.groups;
        },
        curGeneName(){
            return this.$store.state.curGeneName;
        },
        dataset(){
            return this.curData.paramsObj['dataset'];
        },
        chosenNodes(){
            return this.curData.chosenData;
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

    methods:{
        //utils:
        kde:function(kernel, thresholds, data) {
            return thresholds.map(t => [t, d3.mean(data, d => kernel(t - d))]);
        },
        epanechnikov:function(bandwidth) {
            return x => Math.abs(x /= bandwidth) <= 1 ? 0.75 * (1 - x * x) / bandwidth : 0;
        },
        jitter:function(value,group){
            let Intervals = this.curViolinData[group].kde;
            let closestUpper = 1;
            for(let i = 1;i < this.thresholds.length;i++){
                if(value < this.thresholds[i]){
                    closestUpper = i;
                    break;
                }
            }
            let lowerKde = Intervals[closestUpper-1][1];
            let upperKde = Intervals[closestUpper][1];
            let lowerX = this.wiggle(lowerKde);
            let upperX = this.wiggle(upperKde);
            let X = (lowerX + upperX) / 2.0; //TODO:这里可能可以用B-curve优化一下
            //add jitter
            //TODO: 这里使用区间叠加概率随机的方式添加的jitter，可能可以优化一下(优化为越接近边缘概率越低，并且是一种渐变的方式)
            let Part = Math.random();
            let param1 = 0.6
            let param2 = 0.9
            if(Part < param1){
                return (X - Math.random() * 2 * X) * 0.5;
            }
            else if(Part < param2){
                return (X - Math.random() * 2 * X) * 0.7;
            }
            else{
                return (X - Math.random() * 2 * X);
            }
        },
        //main function:
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
        config(){
            const self = this;
            const svg = d3.select('.violin');
            //配置框架
            this.margin = {
                left:40,
                right:25,
                top:25,
                bot:25,
            };

            //辅助计算标签宽度
            const tempG = svg.append('g')
            const tempText = tempG.selectAll('*').data(this.groups.map(v=>v.name)).join('text').text(d=>d).attr('font-size','20px')
            //根据最大的文本尺寸调整大小
            let maxGroupTextSize = 0;
            tempText.each(function(){
                if(this.getBoundingClientRect().width > maxGroupTextSize){
                    maxGroupTextSize = this.getBoundingClientRect().width
                }
            })
            tempG.remove()
            this.unitWidth = Math.max(maxGroupTextSize+25,150)



            this.height = this.$refs.violinContainer.clientHeight * 0.97
            console.log('violin_height:',this.height)
            this.unitPlotMaxWidth = 60
            this.width = this.unitWidth * this.groups.length + this.margin.left + this.margin.right;

            //配置阈值（即各个细分区间的边界值）
            this.thresholds = []
            let IntervalNum = 200;//绘图区间的数目，代表小提琴图曲线的精细度
            for(let i = 0;i <= IntervalNum;i++){
                this.thresholds.push(this.curGeneRange[0] + (this.curGeneRange[1]-this.curGeneRange[0]) * (1.0 * i / IntervalNum));
            }

            //配置数据
            let orderedCells = JSON.parse(JSON.stringify(this.cellData));//注意，这里是一次深copy
            orderedCells.sort(function(cell1,cell2){
                if(self.curGeneExpression[cell1.id] < self.curGeneExpression[cell2.id])
                    return -1;
                else{
                    return 1;
                }
            });
            let dataByGroup = [];
            for(let i = 0;i < this.groups.length;i++){
                dataByGroup.push([
                    this.groups[i].id,[],[]
                ]);
            }
            for(let i = 0;i < orderedCells.length;i++){
                let index = -1;
                for(let j = 0;j < dataByGroup.length;j++){
                    if(dataByGroup[j][0] == orderedCells[i].group){
                        index = j;
                        break;
                    }
                }
                dataByGroup[index][1].push(this.curGeneExpression[orderedCells[i].id]);
                dataByGroup[index][2].push(orderedCells[i].id);
            }
            this.curViolinData= dataByGroup.map(el=>{
                    let q1 = d3.quantile(el[1],0.25);
                    let q3 = d3.quantile(el[1],0.75);
                    let iqr = q3 - q1;
                    let cellInfo = []
                    for(let i = 0;i < el[1].length;i++){
                        cellInfo.push({'id':el[2][i],'expression':el[1][i]});
                    }
                    return ({
                        group:el[0],
                        kde:this.kde(this.epanechnikov(0.2),this.thresholds,el[1]),
                        median:d3.quantile(el[1],0.5),
                        quartile:[d3.quantile(el[1], 0.25),  d3.quantile(el[1], 0.75)],
                        cellInfo:cellInfo,
                        max: Math.max(...cellInfo.map(v=>v.expression)),
                        min: Math.min(...cellInfo.map(v=>v.expression)),
                    });
                }
            );
            // //截断kde中在界限外的部分
            // for(let i = 0;i < this.curViolinData.length;i++){
            //     let newkde = [];
            //     for(let v of this.curViolinData[i].kde){
            //         if(v[0] <= this.curViolinData[i].max && v[0] >= this.curViolinData[i].min){
            //             newkde.push(v)
            //         }
            //     }
            //     this.curViolinData[i].kde = newkde
            // }

            //比例尺配置

            //最大核密度
            let maxKDE = 0;
            for(let i = 0;i < this.curViolinData.length;i++){
                for(let j = 0;j < this.curViolinData[i].kde.length;j++){
                    maxKDE = Math.max(this.curViolinData[i].kde[j][1],maxKDE)
                }
            }

            this.wiggle = d3.scaleLinear()
                            .domain([0, maxKDE * 0.5])
                            .range([0, this.unitPlotMaxWidth * 0.5]);
            this.xScale = d3.scaleBand()
                            .domain(this.curViolinData.map(d => this.groups.find(function(item){
                                return item.id == d.group
                            }).id))
                            .range([this.margin.left, this.width - this.margin.left - this.margin.right]);
            this.yScale = d3.scaleLinear()
                            .domain([this.curGeneRange[0], this.curGeneRange[1]]) //对基因的范围自动适应
                            .range([this.height - this.margin.bot, this.margin.top]);  

            //绘图函数配置
            this.xAxis  = g => {
                g.attr("transform", `translate(0,${this.height - this.margin.bot})`)
                 .call(d3.axisBottom(this.xScale));
            }
            this.yAxis = g => {
                g.attr("transform", `translate(${this.margin.left},0)`)
                 .call(d3.axisLeft(this.yScale))    
            }   
            this.areal = d3
                        .area()
                        .x0(0)
                        .x1(d => this.wiggle(d[1]))
                        .y(d => this.yScale(d[0]))
                        .curve(d3['curveBasis']);
            this.arear = d3
                        .area()
                        .x0(0)
                        .x1(d => -this.wiggle(d[1]))
                        .y(d => this.yScale(d[0]))
                        .curve(d3['curveBasis']);


        },
        drawViolinLegend(){
            const svg = d3.select('.violin');
            svg.style("width",`${this.width}px`);
            svg.style("height",`${this.height}px`);
            let legend = svg.select('.violin-legend')
            if(legend.empty()){
                legend = svg.append('g').classed('violin-legend',true)
            }
            legend.selectAll("*").remove();
            legend.attr("stroke-linejoin","round").attr("stroke-linecap","round");
            legend.style("width",`${this.width}px`);
            legend.style("height",`${this.height}px`);
            legend.append("g").classed('xAxis',true).call(this.xAxis);
            svg.select('.xAxis')
               .selectAll('.tick')
               .select('text')
               .text(d=>{
                return this.groups.find(group=>group.id==d).name
            })
            legend.append("g").classed('yAxis',true).call(this.yAxis);
            //调整文本大小
            legend.select('.xAxis').selectAll('.tick').selectAll('text').attr('font-size','20px')
            legend.select('.yAxis').selectAll('.tick').selectAll('text').attr('font-size','20px')


        },
        drawViolinPlot(){
            const self = this;
            const svg = d3.select('.violin');
            svg.style("width",`${this.width}px`);
            svg.style("height",`${this.height}px`);
            let plot = svg.select('.violin-plot')
            if(plot.empty()){
                plot = svg.append('g').classed('violin-plot',true)
            }
            plot.selectAll("*").remove();
            plot.attr("stroke-linejoin","round").attr("stroke-linecap","round");


            let violins = plot.selectAll(".violin")
                .data(this.curViolinData)
                .enter()
                .append("g")
                .attr("transform", d => `translate(${this.xScale(this.groups.find(e=>{
                    return e.id == d.group
                }).id) + this.xScale.bandwidth() / 2}, 0)`)

            //提琴图填充  
            const fill_opacity = 0.6
            violins.append("path")
                .attr("d",d => this.areal(d.kde))
                .attr("fill", (d) => this.groups.find(e=>d.group == e.id).color)
                .attr("opacity",fill_opacity);
            violins.append("path")
                .attr("d",d => this.arear(d.kde))
                .attr("fill", (d) => this.groups.find(e=>d.group == e.id).color)
                .attr("opacity",fill_opacity);
            


            //提琴图轮廓
            violins.append("path")
                .attr("d",d => this.areal.lineX1()(d.kde))
                .attr("fill-opacity","0")
                // .attr("stroke", "#363636")
                .attr("stroke", (d) => this.groups.find(e=>d.group == e.id).color)
                .attr('stroke-width',4)
            violins.append("path")
                .attr("d",d => this.arear.lineX1()(d.kde))
                .attr("fill-opacity","0")
                // .attr("stroke", "#363636")
                .attr("stroke", (d) => this.groups.find(e=>d.group == e.id).color)
                .attr('stroke-width',4)


            //box plot
            if(this.mode == "Box"){
                violins.append("rect")
                    .attr("x", -8)
                    .attr("y", d => this.yScale(d.quartile[1]))
                    .attr("width", 16)
                    .attr("height", d => this.yScale(d.quartile[0]) - this.yScale(d.quartile[1]))
                    .attr("fill", "black")
                    
            }

            //strip plot
            if(this.mode == "Strip"){
                violins.each(function(d,i){
                    d3.select(this)
                        .selectAll('.Strip')
                        .data(d.cellInfo)
                        .enter()
                        .append('circle')
                        .attr('cx',function(x,j){
                            return self.jitter(x.expression,i);
                        }) //jitter
                        .attr('cy',x =>self.yScale(x.expression))
                        .attr('r',1)
                        .attr('opacity',1)
                        .classed('strip-scatter',true)
                        .on('mouseover',function(e,d){
                            //show info
                            self.infoPanel.show()
                            self.infoPanel.setMessageData({'id':d.id,'value':d.expression})
                            self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)
                        })
                        .on('mousemove',function(e,d){
                            //show info
                            self.infoPanel.show()
                            self.infoPanel.setMessageData({'id':d.id,'value':d.expression})
                            self.infoPanel.setPos(e.clientY - 40,e.clientX + 15)

                        })
                        .on('mouseout',function(){
                            //hidden info
                            self.infoPanel.hidden()

                        })
                });
                this.updateChosenStripScatter();
            }



            
            //box的轴线
            if(this.mode == "Box"){
                violins.append("line")
                    .attr("x1", 0)
                    .attr("x2", 0)
                    .attr("y1", d => this.yScale(d.min))
                    .attr("y2", d => this.yScale(d.max))
                    .attr("stroke", "black")
                    .attr('stroke-width',3)
                violins.append("circle")
                    .attr("cx", 0)
                    .attr("cy", d => this.yScale(d.median))
                    .attr("r", 4)
                    .attr("fill", "white")
                    .attr('stroke','black')
                    .attr('stroke-width',2)
            }


        },
        updateChosenStripScatter(){
            const self = this;
            //如果全部没选，那么全画

            if(this.chosenNodes.length == 0)
            {
                d3.selectAll(".strip-scatter").classed('unchosen',false);
                return;
            }
            const stripScatters = d3.selectAll(".strip-scatter")
                .classed('unchosen',function(d,i){
                    if(self.chosenNodes.indexOf(d.id) != -1){
                        return false;
                    }
                    return true
                })
        },
        async reDraw(){
            eventBus.$emit('GeneExpressionViewRefreshingStart');
            if(this.curData === undefined || this.curData === null)
                return ;
            if(this.curData.cellData === undefined || this.curData.cellData === null){
                return ;
            }
            if(this.curData.cellData.length == 0){//如果数据量为0
                d3.select('.violin').selectAll('*').remove();
            }
            else{                                                                                                                                                                                                                                                                                                                                                                                          
                await this.updateCurGeneInfo();
                this.config();
                this.drawViolinLegend();
                this.drawViolinPlot();
            }
            this.$refs['violin-scroll'].update();
            eventBus.$emit('GeneExpressionViewRefreshingClose')
        },


        menuMounted(_this,root,parent) {
        
        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            // //png保存
            // saveSvgAsPng(this.$refs.violin, "violin.png");
            
            
            //svg保存
            const svgDOM = this.$refs['violin'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "Distribution of Gene Expression View.svg";
            a.click();
            URL.revokeObjectURL(url)

        }
    }, 
    watch:{
        async cellData(){
            if (this.cellData === undefined || this.cellData === null) return;
            this.reDraw();
        },
        curGeneName:{
            deep:true,
            async handler(){
                this.reDraw();
            }
        },
        groups:{
            //主要是监控组名被修改
            deep: true,
            handler(newValue,oldValue) {
                if (oldValue === "undefined" || oldValue === "null") return;
                this.reDraw();
            },
        },
        mode(){
            this.drawViolinPlot();
        },
        chosenNodes(){
            this.updateChosenStripScatter();
        },
        'repaintTag.Violin':{
            handler(){
                this.reDraw();
            }
        }
        
    },
    mounted(){
        this.reDraw()
    }
};
</script>

<style scoped lang="less">
.violin-container{
    height: 100%;
    width:100%;
    
    .violin-plot-container{
        width:100%;
        height:100%;
        
        position: relative;
        /deep/ .el-scrollbar__wrap{
            overflow: hidden;
            position: relative;
        }
        /deep/ .is-horizontal{
            height: 10px;
            .el-scrollbar__thumb{
                background-color:rgb(150, 150, 150);
            }
        }
        /deep/ .is-vertical{
            width: 10px;
            display: none;
            .el-scrollbar__thumb{
            background-color:rgb(150, 150, 150);
            }
        }
        .violin{
            padding: 3px;
        }
    }

    /deep/ .unchosen{
        fill:black;
        opacity: 0;
    }

}



</style>
