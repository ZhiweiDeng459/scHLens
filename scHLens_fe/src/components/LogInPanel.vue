<template>
        <div class="LogInContainer">
            <div style="position:absolute;right:20px;top:10px">
                <el-button  type="text"  @click="$router.push('/tutorial/help')"><b style="font-size:20px">Help</b></el-button>
            </div>
            <b class="LogInTitle">scHLens</b>
            <el-input 
                v-model="loadJobId"
                placeholder="Please input job id if you want to load a exist job"
                style="
                    width:350px;
                    margin-bottom:30px;
                    margin-top:30px;
                "></el-input>
            <div
                style="
                    width:550px;
                    display:flex;
                    justify-content:space-between
                ">
                <el-button
                    type="success"
                    style="width:170px"
                    @click="_createNewJob()"
                    round>Create new Job</el-button>
                <el-button
                    type="primary"
                    style="width:170px;margin-right:7px"
                    @click="_loadExistJob()"
                    round>Load Exist Job</el-button>
                <el-popconfirm
                    confirm-button-text='Yes'
                    cancel-button-text='No'
                    icon="el-icon-info"
                    icon-color="orange"
                    confirm-button-type="Primary"
                    cancel-button-type="Primary"
                    @confirm="handleClickJobLoaderButton"
                    title="Warning! Uploading the job will overwrite the existing job with the same ID. Do you still want to proceed?">
                    <el-button slot="reference" type="primary" style="width:170px" round>Import Saved Job</el-button>
                </el-popconfirm>

            </div>

            <div
                style="width: 100%;
                display: flex;
                justify-content: flex-end;">


                <el-popover
                    placement="bottom"
                    trigger="hover">
                    <div>

                        <el-table 
                            :data="SampleJobInfo" 
                            :show-header="false" 
                            @row-click="handleClickJobRow"
                            :max-height="550"
                            border>
                            <el-table-column width="300" property="name" label="name"></el-table-column>
                        </el-table>

                    </div>
                    <el-button slot="reference" type="text" style="margin-right: 0px;margin-top:10px;"><b style="font-size: 18px;">>>> Use Examples?</b></el-button>
                </el-popover>

                <input type="file" ref="JobLoader" @change="handleJobLoaderChangeFile" style="display: none" /> <!--真实的job上传器-->

            </div>

        </div>

</template>

<script>
import axios from "axios";
import {createNewJob,loadExistJob,createNewSocket,fetchSampleJobInfo,loadSampleJob,uploadJob} from "@/utils/interface"
import eventBus from '@/utils/eventBus.js'
import { Form, Tabs, TabPane, Loading, Drawer,Message,Popconfirm} from "element-ui";



export default {
    name: "LogInPanel",

    props:['closeLogIn','setJobId'],

    data(){
        return {
            'loadJobId':'',
            'SampleJobInfo':[],//Sample Job的信息 {'name':'3k-PBMCs','id':'3k-PBMCs'}
        }
    },

    computed:{
        JobId(){
            return this.$store.state.JobId;
        },
        socketIns(){
            return this.$store.state.socketIns;
        }
    },

    methods:{
        _createNewJob(){
            const loading = Loading.service({ fullscreen: true });
            createNewJob()
                .then((response) => {
                    let newJobId = response.data
                    this.$store.commit("setJobId",newJobId) //更新JobId
                    this.$store.commit("updateDatasets",newJobId) //更新数据集信息
                    //创建与后台的socket连接
                    this.$store.commit("setSocket", createNewSocket(newJobId));
                    //后台消息管理
                    this.socketIns.on('general_info',(data)=>{//监听后台消息，并呈现
                        this.$notify({
                            title: data['title'],
                            message: data['content'],
                            type:data['type'],
                            duration:15000,
                        })
                    })



                    //关闭对话框
                    this.closeLogIn()
                    //关闭loading
                    loading.close();

                    
                })
                .catch((err)=>{
                    //关闭loading
                    loading.close();
                    this.$message({
                        'message':'Unable to connect to the server',
                        'type':'error',
                        'showClose':true,
                    })
                    console.log(err)
                })

        },
        _loadExistJob(){//外部函数
            const loading = Loading.service({ fullscreen: true ,text:'It may take several minutes or even tens of minutes, please wait...'});
            if(this.loadJobId == ''){
                loading.close();
                this.$message({
                    'message':'Unexist Job ID',
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            loadExistJob(this.loadJobId)//网络接口
                .then((response)=>{
                    console.log('load res',response)
                    if(response.data == 'unexist'){ //输入了不存在的JobId
                        loading.close();
                        this.$message({
                            'message':'Unexist Job ID',
                            'type':'error',
                            'showClose':true,
                        })       
                    }
                    else{  //存在的Job
                        if(Object.keys(response.data).length != 0)
                            eventBus.$emit("receivedTreeData",response.data)
                        this.$store.commit("setJobId",this.loadJobId)
                        this.$store.commit("updateDatasets",this.loadJobId)
                        //创建socket连接
                        this.$store.commit("setSocket", createNewSocket(this.loadJobId));
                        //关闭对话框
                        this.closeLogIn()
                        loading.close();
                    }

                })
                .catch((err)=>{
                    console.log(err);
                    this.$message({
                        'message':'Unable to connect to the server',
                        'type':'error',
                        'showClose':true,
                    })
                    loading.close();
                })

        },
        getSampleJobInfo(){
            fetchSampleJobInfo().then((res)=>{
                this.SampleJobInfo = res.data.slice().reverse();
            }).then((err)=>{
                console.log(err)
            })
        },
        handleClickJobRow(row,column,event){
            let sampleId = row.id;
            const loading = Loading.service({ 
                fullscreen: true, 
                text: 'It may take several minutes or even tens of minutes, please wait...',
            });
            loadSampleJob(sampleId).then((res)=>{
                loading.close()
                let newJobId = res.data['newJobId']
                this.loadJobId = newJobId
                this._loadExistJob();

            }).catch((err)=>{
                loading.close()
                this.loadJobId = ''
                this.$message({
                    'message':'Load failed!',
                    'type':'error',
                    'showClose':true,

                })
                return;
            })
        },
        handleClickJobLoaderButton(){
            this.$refs['JobLoader'].click();
        },
        handleJobLoaderChangeFile(event){
            const loading = Loading.service({ fullscreen: true ,text:'Importing... This may take some time.'});
            let form = new FormData()
            form.append('file',event.target.files[0])
            uploadJob(form).then(res=>{
                this.$message({
                    'message':'UPLOAD COMPLETE',
                    'type':'success',
                    'showClose':true,
                })
                loading.close();
                //加载最新的job
                let uploadJobId = res.data['uploadJobId']
                this.loadJobId = uploadJobId
                this._loadExistJob();
            },err=>{
                this.$message({
                    'message':'UPLOAD FAILED',
                    'type':'error',
                    'showClose':true,
                })
                loading.close();
            })
        }
    },
    mounted(){
        //获取Sample Job信息
        this.getSampleJobInfo();
    }
}
</script>

<style scoped lang="less">
    .LogInContainer{
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 550px;
        height: 250px;
        // background-color: #24292f;
        .LogInTitle{
            font-size: 50px;
            font-family: AslinaBold;
            color: #90e36b;
            margin-top:15px;
            margin-bottom:20px;  
        }
    }
    /deep/ .el-table__row{
        cursor: pointer;
    }

</style>