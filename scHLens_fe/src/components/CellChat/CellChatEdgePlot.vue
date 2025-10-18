<template>
  <div>
    <SelfContextMenu
            :items = menuItems
            :_mounted = menuMounted
        />
    <svg ref="CellChatEdgePlot" id="CellChatEdgePlot" style="background-color: white">
    </svg>
  </div>
</template>

<script>
import * as d3 from "d3";
import {saveSvgAsPng} from 'save-svg-png-ext'
import SelfContextMenu from "@/components/SelfContextMenu"
import eventBus from "@/utils/eventBus.js"

export default {
    name:'CellChatEdgePlot',
    props:['setInteractionTable','dataMode'],
    components:{
        'SelfContextMenu':SelfContextMenu
    },

    data() {
        return {
            mode:'count',
            menuItems:[
                {
                    'name':'Save this Image',
                    'icon':'icons/save_as_image.svg',
                    'callback':()=>{
                        this.saveToFile();
                    }
                },
                {
                    'name':'Save interaction table',
                    'icon':'icons/download.svg',
                    'callback':()=>{
                        const self = this
                        let interaction_map = {}
                        let groupids = self.groups.map(item=>item['id']) 
                        for(let i = 0;i < groupids.length;i++){
                            interaction_map[groupids[i]] = {} 
                            for(let j = 0;j < groupids.length;j++){
                                interaction_map[groupids[i]][groupids[j]] = []
                            }
                        }
                        for(let i = 0;i < self.CC['source'].length;i++){
                            let source = self.CC['source'][i]
                            let target = self.CC['target'][i]
                            let ligand = self.CC['ligand'][i]
                            let receptor = self.CC['receptor'][i]
                            let score = self.CC['prob'][i]
                            interaction_map[source][target].push({
                                'source':self.groups.find(v=>v.id==source).name,
                                'target':self.groups.find(v=>v.id==target).name,
                                'ligand':ligand,
                                'receptor':receptor,
                                'score':score.toPrecision(4)
                            })
                        }

                        let table_data = []
                        for(let i = 0;i < groupids.length;i++){
                            for(let j = 0;j < groupids.length;j++){
                                table_data.push(...interaction_map[groupids[i]][groupids[j]])
                            }
                        }

                        if(table_data.length==0){
                            this.$message({
                                message: 'No interaction data to export!',
                                type: 'error'
                            });
                            return;
                        }

                        // 1. 将对象数组转换为 CSV 字符串
                        const headers = Object.keys(table_data[0]).join(",") + "\n";
                        const rows = table_data.map(obj => Object.values(obj).join(",")).join("\n");
                        const csvContent = headers + rows;

                        // 2. 创建一个 Blob 对象
                        const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });

                        // 3. 创建下载链接并触发下载
                        const link = document.createElement("a");
                        link.href = URL.createObjectURL(blob);
                        link.download = "interaction.csv";
                        document.body.appendChild(link);
                        link.click();
                        document.body.removeChild(link);

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
        dataMode(){
            this.reDraw()
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

            if(this.CC===undefined || this.CC===null) return;
            if(Object.keys(this.CC).length===0) return;
            //添加zoom层
            const zoomSVG = svg.append('g')
                               .attr("class", "cellchat-edge-zoomLayer")
            //定义zoom行为
            const zoom = d3.zoom()
                               .scaleExtent([0.2, 20])  // 设置缩放范围，最小 0.5 倍，最大 5 倍
                               .on("zoom", (event) => {
                                    zoomSVG.attr("transform", event.transform);  // 平移缩放
                                });
            svg.call(zoom)


            //计算边界
            let cellCount = 0; //细胞总数
            this.groups.forEach(element => {
                cellCount += element['size']
            });
            
            const maxGroupSize = Math.max(...this.groups.map(v=>v['size'])) //最大的group的细胞数
            const minGroupSize = Math.min(...this.groups.map(v=>v['size'])) //最小的group的细胞数
 
            const maxDistance = 0.5 * size / (1 + Math.sin(Math.PI / this.groups.length)); //保证聚类圆不相交，不溢出的又可理论接受的聚类圆心到点的距离，TODO 这里设置为默认距离

            const maxRadius = 0.5 *  maxDistance * Math.sin(Math.PI / this.groups.length); //在maxDistance的前提下，理论可接受的最大聚类半径。如果设定距离小于maxDistance，那么半径也要适当缩短

            const RadiusRange = [maxRadius * 0.2,maxRadius * 0.7]; //聚类圆的上下限，最大size对应上限，最小size对应下限
            const RadiusRangeScaler = d3.scaleLinear()
                                        .domain([minGroupSize,maxGroupSize])
                                        .range(RadiusRange)
            

            //整理数据
            let rawData = this.CC[this.dataMode];
            let names =  this.groups.map(item=>item['id'])
            let attachGroupData = this.groups.map((item,i)=>{ //groups的辅助信息
                return {
                    'id':item['id'],
                    'angle': i * 2 * Math.PI / self.groups.length,//圆心所在的角度
                    'x' : centerX + Math.cos( i * 2 * Math.PI / self.groups.length ) * maxDistance,//圆心的横坐标
                    'y' : centerY + Math.sin( i * 2 * Math.PI / self.groups.length ) * maxDistance, //圆心的纵坐标
                    'r' : RadiusRangeScaler(this.groups[i].size), //圆的半径,按照
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

            const maxValue = Math.max(...drawData.map(v=>v.value))
            const minValue = Math.min(...drawData.map(v=>v.value))

            const StrokeRange = [maxRadius * 0.08, maxRadius * 0.6]
            const StrokeRangeScaler = d3.scaleLinear()
                                        .domain([minValue,maxValue])
                                        .range(StrokeRange)

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

            //绘制圆
            const node = zoomSVG.append('g');
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

            const group_text = zoomSVG.append('g')


            //定义箭头
            zoomSVG.append('defs')
                .append('marker')
                .attr('id','cellchat-edge-arrowhead')
                .attr('viewBox','0 0 10 10')
                .attr('refX',10)
                .attr('refY',5)
                .attr('orient','auto')
                .attr('markerWidth',0.4 * maxRadius)
                .attr('markerHeight',0.4 * maxRadius)
                .attr('markerUnits','userSpaceOnUse')
                .append('path')
                .attr('d','M 0 0 L 10 5 L 0 10 z')
                .attr('fill','black')
                .style('stroke','none');

            //绘制连接线
            const link = zoomSVG.append('g')
            link.selectAll('path')
                .data(drawData)
                .join('path')
                .classed('cell-chat-edge-link',true)
                .attr('d',function(d,i){
                    //计算端点的位置
                    let source = attachGroupData.find(v=>v.id == d.source)
                    let target = attachGroupData.find(v=>v.id == d.target)
                    let distance = Math.sqrt( 
                        Math.pow(source.x - target.x,2) +  Math.pow(source.y - target.y,2)
                        )
                        
                    
                    const link = d3.line()
                                .curve(d3.curveBasis)
                                .x(d=>d.x)
                                .y(d=>d.y)

                    if(d.source != d.target){//起点终点不同
                        let diver_dis = maxDistance * 0.1; //向外侧偏移
                        let _target = { //位于两个圆心连接线上，在source圆边上的点
                            'x':target.x + (1.0 * target.r / distance) * (source.x - target.x),
                            'y':target.y + (1.0 * target.r / distance) * (source.y - target.y),
                        }

                        let _center = { //连线中点
                            'x': 0.5 * (source.x + target.x),
                            'y': 0.5 * (source.y + target.y)
                        }

                        let _center_c = { //连接线的中段控制点，用于产生平滑的弯曲
                            'x':_center.x + (source.y - target.y) / distance * diver_dis,
                            'y':_center.y + (target.x - source.x) / distance * diver_dis,
                        }
                        return link([source,_center_c,_target])
                    }
                    else{ //起点终点相同
                        let two_circle_distance = Math.sqrt( //小圆到大圆中心的距离
                            Math.pow(centerX - target.x,2) +  Math.pow(centerY - target.y,2)
                        )
                        let _target = { //线段终点，位于两个圆心连接线上，在target圆边上的点
                            'x':target.x + (1.0 * target.r / two_circle_distance) * (target.x - centerX),
                            'y':target.y + (1.0 * target.r / two_circle_distance) * (target.y - centerY),
                        }
                        let _source = _target; //起点和终点重合

                        let dx = _source.x - centerX;
                        let dy = _source.y - centerY;
                        let len = Math.sqrt(dx*dx + dy*dy);
                        dx /= len; dy /= len;

                        let offset = maxRadius * 1.5; // 外扩距离
                        let rx = maxRadius * 0.8; 
                        let ry = maxRadius * 1.2;


                        // 椭圆上的两个对称控制点
                        let c1 = { x: _source.x + dx*offset - dy*ry, y: _source.y + dy*offset + dx*rx };
                        let c2 = { x: _source.x + dx*offset + dy*ry, y: _source.y + dy*offset - dx*rx };

                        return link([_source,c1,c2,_target]);//用两个控制点构造贝塞尔曲线
                    }

                })
                .attr('stroke',function(d){
                    return self.groups.find(v=>v.id == d.source).color
                })
                .attr('stroke-width',(d)=>{
                    return `${StrokeRangeScaler(d.value)}px`
                })
                .attr('opacity',0.5)
                .attr('fill',function(d){
                    return 'none'
                })
                .attr('display',d=>{
                    if(d.value == 0)
                        return 'none'
                    else
                        return 'block'
                })
                .attr('marker-end',"url(#cellchat-edge-arrowhead)")
                .on("mouseover",function(e,d){
                    //突出显示
                    d3.select(this)
                        .attr("opacity",1)
                    //显示信息
                    self.infoPanel.show()
                    if(self.dataMode=='weight'){
                        self.infoPanel.setMessageData({
                            'Source':self.groups.find(v=>v.id == d.source).name,
                            'Target':self.groups.find(v=>v.id == d.target).name,
                            'Strength':rawData[d.source][d.target].toPrecision(3),
                        })
                    }
                    else{//count
                        self.infoPanel.setMessageData({
                            'Source':self.groups.find(v=>v.id == d.source).name,
                            'Target':self.groups.find(v=>v.id == d.target).name,
                            'Count':rawData[d.source][d.target],
                        })
                    }
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 30)
                    d3.select(this).raise()//放到最顶层
                })
                .on("mousemove",function(e,d){
                    d3.select(this)
                        .attr("opacity",1)
                    //显示信息
                    self.infoPanel.show()
                    if(self.dataMode=='weight'){
                        self.infoPanel.setMessageData({
                            'Source':self.groups.find(v=>v.id == d.source).name,
                            'Target':self.groups.find(v=>v.id == d.target).name,
                            'Strength':rawData[d.source][d.target].toPrecision(3),
                        })
                    }
                    else{//count
                        self.infoPanel.setMessageData({
                            'Source':self.groups.find(v=>v.id == d.source).name,
                            'Target':self.groups.find(v=>v.id == d.target).name,
                            'Count':rawData[d.source][d.target],
                        })
                    }
                    self.infoPanel.setPos(e.clientY - 40,e.clientX + 30)

                })
                .on("mouseout",function(e,d){
                    d3.select(this)
                        .attr("opacity",0.5)
                    self.infoPanel.hidden()

                })
                .on("click",function(e,d){
                    let table_data = interaction_map[d.source][d.target]

                    self.setInteractionTable(table_data)
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

        
            
            //绘制文本
            // group_text.selectAll('text')
            //     .data(this.groups)
            //     .join('text')
            //     .attr('x',function(d,i){
            //         return attachGroupData[i].x;
            //     })
            //     .attr('y',function(d,i){
            //         if(attachGroupData[i].y < centerY)
            //             return attachGroupData[i].y - attachGroupData[i].r - 10;
            //         else{
            //             return attachGroupData[i].y + attachGroupData[i].r + 10;
            //         }
                    
            //     })
            //     .attr('text-anchor','middle')
            //     .attr('dominant-baseline','middle')
            //     .attr('font-size','14px')
            //     .style('font-family','YaHei')
            //     .style('font-weight','bold')
            //     .text(d=>d.name)
            group_text.selectAll('text')
                .data(this.groups)
                .join('text')
                .attr('x', function(d, i) {
                    const data = attachGroupData[i];
                    // Calculate the vector from center (centerX, centerY) to the circle center (data.x, data.y)
                    const dx = data.x - centerX;
                    const dy = data.y - centerY;
                    const dist = Math.sqrt(dx * dx + dy * dy);

                    // Normalize the vector and scale it by the circle's radius (data.r)
                    // This calculates the x-coordinate of the point on the circle's edge
                    // along the line from the center.
                    // We add a small buffer (e.g., 5) to move the text slightly outside.
                    const buffer = 5;
                    return data.x + (dx / dist) * (data.r + buffer);
                })
                .attr('y', function(d, i) {
                    const data = attachGroupData[i];
                    // Calculate the vector from center (centerX, centerY) to the circle center (data.x, data.y)
                    const dx = data.x - centerX;
                    const dy = data.y - centerY;
                    const dist = Math.sqrt(dx * dx + dy * dy);

                    // Normalize the vector and scale it by the circle's radius (data.r)
                    // This calculates the y-coordinate of the point on the circle's edge
                    // along the line from the center.
                    const buffer = 5;
                    return data.y + (dy / dist) * (data.r + buffer);
                })
                .attr('text-anchor', 'start') // Anchor the text to the *start* of the label
                .attr('dominant-baseline', 'middle') // Keep it vertically centered on the line
                .attr('font-size', '14px')
                .style('font-family', 'YaHei')
                .style('font-weight', 'bold')
                .attr('transform', function(d, i) {
                    const data = attachGroupData[i];
                    const x = data.x;
                    const y = data.y;

                    // Calculate the angle (in degrees) from the center (centerX, centerY) to the circle center (x, y)
                    // atan2(dy, dx) returns radians in the range [-pi, pi]
                    let angle = Math.atan2(y - centerY, x - centerX) * (180 / Math.PI);

                    // Adjust the angle:
                    // The angle calculated is for the line from the center (centerX, centerY) to the circle center (x, y).
                    // We want the text to read *outward*.
                    // If the angle is in the left hemisphere (90 to 270 degrees), we flip the text 180 degrees
                    // and change the text-anchor to 'end' so it reads correctly.

                    if (angle > 90 || angle < -90) {
                        angle += 180; // Flip the text for the left side
                        // For the flipped text, we'll want to adjust the text-anchor:
                        d3.select(this).attr('text-anchor', 'end');
                    } else {
                        // Ensure text-anchor is 'start' for the right side
                        d3.select(this).attr('text-anchor', 'start');
                    }

                    // Apply the rotation transform around the new (x, y) position of the text
                    const textX = d3.select(this).attr('x');
                    const textY = d3.select(this).attr('y');

                    return `rotate(${angle}, ${textX}, ${textY})`;
                })
                .text(d => d.name);

        },
        saveToFile(){
            /**
             * 保存视图为文件
             */
            // //png保存
            // saveSvgAsPng(this.$refs.dotplot, "dotplot.png");
            //svg保存
            const svgDOM = this.$refs['CellChatEdgePlot'];
            const svgData = new XMLSerializer().serializeToString(svgDOM);
            
            const blob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"})
            const url = URL.createObjectURL(blob)

            const a = document.createElement("a")
            a.href = url;
            a.download = "Cell Communication View - EdgePlot.svg";
            a.click();
            URL.revokeObjectURL(url)

            // const svgDOM = this.$refs['CellChatEdgePlot'];
            // const clone = svgDOM.cloneNode(true); // 克隆整个SVG
            // const g = clone.querySelector("g");

            // // 移除缩放变换
            // g.removeAttribute("transform");

            // // 导出克隆版
            // const svgData = new XMLSerializer().serializeToString(clone);
            // const blob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" });
            // const url = URL.createObjectURL(blob);
            // const a = document.createElement("a");
            // a.href = url;
            // a.download = "Cell Communication View - EdgePlot.svg";
            // a.click();
            // URL.revokeObjectURL(url);

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
#CellChatEdgePlot{
    height:100%;
    width:100%;
    cursor:move;
}

.cell-chat-edge-link{
    cursor: pointer;
}

</style>