<!--
  bug:
    1.无法获取大小
    2.如果不实现调用setMessageData，那么位置无法更新
-->
<template>
  <div ref="info-panel-container" v-show="showFlag" class="info-panel-container" :style=InfoPanelStyle>
    
    <!--Gene Table-->
    <el-table
      :data="tableData"
      height="300"
      @row-click="handleTableRowClick"
      border>
      <el-table-column
        prop="gene_name"
        :label="gene_title">
        <!-- eslint-disable vue/no-unused-vars -->
        <template slot="header" slot-scope="unused">
          <div
            style="
              width: 100%;
              height: 20px;
              display: flex;
              justify-content: space-between;
              align-items: center;
              cursor:move"
              @mousedown="startDrag">
            <a style="font-size: 16px;color:black">{{gene_title}}</a>
            <el-button circle type="danger" @click="handlePanelCloseButton" icon="el-icon-close" style="width:16px;height:16px; padding: 0px;display: flex;align-items: center;justify-content: center;"></el-button>
          </div>
        </template>
      </el-table-column>
    </el-table>
    

  </div>
</template> 



<script>

import {Table,TableColumn} from 'element-ui'
import Vue from 'vue'

Vue.component(Table.name,Table)
Vue.component(TableColumn.name,TableColumn)

export default {
  name: 'AnnoRecomAddiPanel',
  data () {
    return {

      tableData:[//用于显示的message数据

      ],
      gene_title:'',
      showFlag:false,
      curPos:{
        'x':0,
        'y':0
      },
      InfoPanelStyle:{
        'top':'0px',
        'left':'0px',
      },
      //拖拽相关
      isDragging:false,
      dragOffset:{
        'x':0,
        'y':0
      }
    }
  },

  computed:{
    curGeneName(){
      return this.$store.state.curGeneName
    }
  },

  methods:{
    setGeneData(gene_title,gene_list){
      this.tableData = gene_list.map(gene=>{
        return {
          'gene_name':gene
        }
      })
      this.gene_title = gene_title
    },
    show(){
      this.showFlag = true
    },
    hidden(){
      this.showFlag = false
    },
    setPos(top,left){
      this.curPos.x = left
      this.curPos.y = top
      this.InfoPanelStyle['top'] = `${this.curPos.y}px`
      this.InfoPanelStyle['left'] = `${this.curPos.x}px`
    },
    handlePanelCloseButton(){
      this.hidden()      
    },
    handleTableRowClick(row){//表格行被点击事件
      let gene = row.gene_name
      //自动添加到currnetGene里面去
      if(this.curGeneName.includes(gene)){//基因已经有了该基因
          this.$message({
              'message':'This Gene has been added',
              'type':'error',
              'showClose':true,
          })
          return;
      }
      //上传基因修改事件
      this.$store.commit("addToCurGeneName", gene);

    },
    //拖拽功能
    startDrag(e){
      this.isDragging = true

      this.dragOffset.x = e.clientX - this.curPos.x;
      this.dragOffset.y = e.clientY - this.curPos.y;

      document.addEventListener("mousemove", this.onDrag);
      document.addEventListener("mouseup", this.stopDrag);
    },
    onDrag(e) {
      if (this.isDragging) {

        this.setPos(e.clientY - this.dragOffset.y,e.clientX - this.dragOffset.x)

      }
    },
    stopDrag() {
      this.isDragging = false;
      document.removeEventListener("mousemove", this.onDrag);
      document.removeEventListener("mouseup", this.stopDrag);
    },

  },
  mounted(){
    document.addEventListener('click',(e)=>{ //任意的单击事件都会关闭菜单栏
            this.hidden()
    })
  }

}
</script>



<style scoped>
  .info-panel-container{
    position:absolute;
    z-index:9999;
    width:220px;
    /* 边框加深 */
    /* border: 1px solid #555;  */

    /* 阴影和悬浮感 */
    box-shadow: 0 8px 16px rgba(0, 0, 0, 0.2), /* 主阴影，向下和右偏移，模糊更大，颜色稍深 */
                0 3px 6px rgba(0, 0, 0, 0.1);   /* 次级阴影，偏移较小，模糊较小，颜色稍浅，增加层次感 */

  }

</style>