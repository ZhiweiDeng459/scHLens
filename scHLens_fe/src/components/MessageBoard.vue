
<template>


  <el-dialog
    :visible.sync="visFlag"
    width="800px"
    :show-close="true"
    :lock-scroll="true"
    :modal="true"
    :close-on-click-modal="false">

    <b style="font-size:23px" slot="title">{{title}}</b>

    <div ref="message-board-container" class="message-board-container" :style=InfoPanelStyle>

      <div style="width:100%">
        <div v-for="(m,index) in messageData" :key="index" style="display:flex;align-items:center;margin-bottom: 10px;">
          <a style="font-size:18px;color:black"><b style="font-size:21px;color:black">{{m[0]}}：</b>{{m[1]}}</a>
        </div>
      </div>

    </div>

  </el-dialog>

</template>



<script>

export default {
  name: 'MessageBoard',
  data () {
    return {

      title:'',

      /**
       * messageData:
       * [
       *    ['attr1',attr1],
       *    ['attr2',attr2],
       *    ...
       * ]
       * 
       */
      messageData:[//用于显示的message数据

      ],

      visFlag:false,
      MessageBoardStyle:{

      },

      InfoPanelStyle:{

      },
    }
  },

  methods:{
    setTitle(title){
      this.title = title
    },

    setMessageData(message){
      // message = {
      //   'a':'asa',
      //   'b':'vvv'
      // }
      this.messageData.length = 0;
      for(let key in message){
        this.messageData.push([key,message[key]])
      }
    },
    show(){
      this.visFlag = true
    },
    hidden(){
      this.visFlag = false
    },

  }

}
</script>



<style scoped lang="less">
.message-board-container{
  padding:10px 20px;
  border-radius:4px;
  border:2px solid gray;
  background-color: white;
  overflow-y: auto;
  height:400px;
}

/deep/ .el-dialog__body{//流水线对话框的主体（除去标题）
  background-color:#F5F5F5;
  border-top: 2px solid lightgray;

}

</style>