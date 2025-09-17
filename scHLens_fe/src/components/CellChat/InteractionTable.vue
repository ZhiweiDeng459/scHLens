<template>
  <div>
    
    <!--Interaction Table-->
    <el-table
      :data="tableData"
      @cell-click="handleTableCellClick"
      size="mini"
      height="240px"
      style="width:250px"
      class="i-table"
      
      border>
      <el-table-column
        prop="ligand"
        label="Ligand"
        width="80px"
        header-align="center"
        class-name="itable-ligand-columm"
        align="center">
      </el-table-column>
      <el-table-column
        prop="receptor"
        label="Receptor"
        width="80px"
        header-align="center"
        class-name="itable-receptor-columm"
        align="center">
      </el-table-column>
      <el-table-column
        prop="score"
        label="Score"
        header-align="center"
        width="80px"
        align="center">
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
  name: 'InteractionTable',
  data () {
    return {

      tableData:[//用于Interaction数据

      ],
    }
  },

  computed:{
    curGeneName(){
      return this.$store.state.curGeneName
    }
  },

  methods:{
    setData(interactionData){
      this.tableData = interactionData
    },


    handleTableCellClick(row, column, cell, event){//表格行被点击事件
      
      if(column.property == 'ligand'){
        let genelist = row.ligand.split('_')
        //统计多少个基因已经在curGene中
        let count = 0
        for(let gene of genelist){
          if(this.curGeneName.includes(gene)){
            count += 1
          }
        }
        if(count == genelist.length && count==1){//说明所有基因都已经在curGene中
          this.$message({
              'message':'This gene has been added before',
              'type':'error',
              'showClose':true,
          })
        }
        else if(count == genelist.length && count > 1){//说明所有基因都已经在curGene中
          this.$message({
              'message':'These genes have been added before',
              'type':'error',
              'showClose':true,
          })
        }
        else if(count == 1){//说明一个基因已经在curGene中
          this.$message({
              'message':`1 gene has been added before`,
              'type':'warning',
              'showClose':true,
          })
        }
        else if(count > 1){//说明部分基因已经在curGene中
          this.$message({
              'message':`${count} genes have been added before`,
              'type':'warning',
              'showClose':true,
          })
        }
        //提交基因
        this.$store.commit("addGroupToCurGeneName", genelist);

      }
      else if(column.property == 'receptor'){
        let genelist = row.receptor.split('_')
        //统计多少个基因已经在curGene中
        let count = 0
        for(let gene of genelist){
          if(this.curGeneName.includes(gene)){
            count += 1
          }
        }
        if(count == genelist.length && count==1){//说明所有基因都已经在curGene中
          this.$message({
              'message':'This gene has been added before',
              'type':'error',
              'showClose':true,
          })
        }
        else if(count == genelist.length && count > 1){//说明所有基因都已经在curGene中
          this.$message({
              'message':'These genes have been added before',
              'type':'error',
              'showClose':true,
          })
        }
        else if(count == 1){//说明一个基因已经在curGene中
          this.$message({
              'message':`1 gene has been added before`,
              'type':'warning',
              'showClose':true,
          })
        }
        else if(count > 1){//说明部分基因已经在curGene中
          this.$message({
              'message':`${count} genes have been added before`,
              'type':'warning',
              'showClose':true,
          })
        }
        //提交基因
        this.$store.commit("addGroupToCurGeneName", genelist);
      }

    },


  },
  mounted(){

  }

}
</script>



<style scoped>
  .i-table{
    /* 边框加深 */
    border: 1px solid lightgray; 

  }
  /deep/ .itable-ligand-columm{
    cursor: pointer;
  }
  /deep/ .itable-receptor-columm{
    cursor: pointer;
  }
</style>