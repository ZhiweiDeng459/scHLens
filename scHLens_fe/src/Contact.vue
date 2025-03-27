<template>
    <div class="contact-container">
        <div class="content-box AboutUs">
            <div class="content-header">
                <div style="flex:1 1 0">
                    <a class="content-title">About Us</a>
                </div>   
            </div>
            <div class="about-unit">
                <div class="about-unit-header">
                    <img src="icons/location.svg" class="about-unit-header-icon"/>
                    <b class="about-unit-header-text">Location</b>
                </div>
                <ul>
                    <li style="margin-bottom: 15px;"><a>Information building #429A, Central South University, Yuelu District, Changsha, Hunan, China</a></li>
                    <li style="margin-bottom: 15px;"><a>Computer building #412, Central South University, Yuelu District, Changsha, Hunan, China</a></li>

                </ul>
            </div>


            <hr width="90%" size="2" color="lightgray" >

            <div class="about-unit">
                <div class="about-unit-header">
                    <img src="icons/email.svg" class="about-unit-header-icon"/>
                    <b class="about-unit-header-text">Email</b>
                </div>
                <ul style="margin-bottom: 15px;">
                    <li style="margin-bottom: 15px;"><a>rqzheng@csu.edu.cn</a></li>
                    <li style="margin-bottom: 15px;"><a>zhiweideng@csu.edu.cn</a></li>

                </ul> 
            </div>

        </div>
        <div class="content-box SendBox">
            <div class="content-header">
                <div style="flex:1 1 0">
                    <a class="content-title">Leave us a message!</a>
                </div>
            </div>
            <div class="send-box-intro">
                <a>If you have any questions/suggestions or if you want to integrate your methods into scHLens, please fill the form below, we will reply you as soon as possible. We look forward to hearing from you!</a>
            </div>
            <el-form v-model="sendForm" label-position="left" label-width="90px" class="send-box-form">
                <el-form-item label="First Name">
                    <el-input size="mini" v-model="sendForm.first_name" placeholder="First Name"></el-input>
                </el-form-item>
                <el-form-item label="Last Name">
                    <el-input size="mini" v-model="sendForm.last_name" placeholder="Last Name"></el-input>
                </el-form-item>
                <el-form-item label="Email">
                    <el-input size="mini" v-model="sendForm.email" placeholder="Email"></el-input>
                </el-form-item>
                <el-form-item label="Message">
                    <el-input size="mini" type="textarea" :rows="7" v-model="sendForm.message"></el-input>
                </el-form-item>
                <el-form-item>
                    <el-button type="primary" @click="submitMessageForm">Send</el-button>
                </el-form-item>

            </el-form>
        </div>
    </div>
</template>

<script>
import Vue from "vue";
import { Form,FormItem,Loading} from "element-ui";
import {sendMessage} from "@/utils/interface"
Vue.component(Form.name, Form);
Vue.component(FormItem.name, FormItem);

export default {
    name: "Contact",
    components: {

    },
    data(){
        return {
            sendForm:{
                first_name:'',
                last_name:'',
                email:'',
                message:'',
            }
        }
    },
    computed: {

    },
    methods:{
        submitMessageForm(){
            
            const loading = Loading.service({ fullscreen: true });
  
            let sendForm = this.sendForm

            //上传校验
            if(sendForm.first_name == '' || sendForm.last_name==''){
                this.$message({
                    'message':'Please enter a valid first name and last name',
                    'type':'error',
                    'showClose':true
                })
                loading.close()
                return;
            }
            if(sendForm.email == ''){
                this.$message({
                    'message':'Please enter a valid email address',
                    'type':'error',
                    'showClose':true,
                })
                loading.close()
                return;
            }



            sendMessage(sendForm)
                .then((res)=>{
                    this.$message({
                        'message':'Your message has been sent successfully',
                        'type':'success',
                        'showClose':true,
                    })
                    this.sendForm.first_name = ''
                    this.sendForm.last_name = ''
                    this.sendForm.email = ''
                    this.sendForm.message = ''
                    loading.close()
                })
                .catch((err)=>{
                    this.$message({
                        'message':'Runtime error',
                        'type':'error',
                        'showClose':true,
                    })
                    console.log(err)
                    loading.close
                })
        }
    },
    mounted(){

    }
};
</script>



<style scoped lang="less">
    .contact-container{
        background-color:#F5F5F5;
        display: flex;
        flex: 1 1 0;
        flex-direction: row;
        align-items: center;
        justify-content: center;
        overflow-x:hidden;
        overflow-y:auto;
        width: 100%;
        .content-box{
            background-color: white;
            margin: 10px;
            border: 2px solid rgb(200, 200, 200);
            border-radius: 4px;
            box-shadow: 0px 0px 9px lightgray;
        }
        a{
            color:#434343
        }
        .AboutUs{
            width:380px;
            height:600px;
            display: flex;
            flex-direction: column;
            .about-unit{
                padding:15px 15px;
                .about-unit-header{
                    margin-bottom:20px;
                    .about-unit-header-icon{
                        width:20px;
                        height:20px;
                        margin-right:10px;
                    }
                    .about-unit-header-text{
                        font-size:18px;
                    }
                }
            }
        }
        .SendBox{
            width:900px;
            height: 600px;
            display: flex;
            flex-direction: column;
            .send-box-intro{
                padding:30px 30px;
            }
            .send-box-form{
                padding:10px 30px;
            }
        }
        .content-header{
            padding-left:15px;
            flex:0 0 40px;
            display: flex;
            align-items: center;
            background-color: #24292f;
            border:2px solid #24292f;
            .content-title{
                // font-family:YaHei;
                color:white;
                font-size:20px;
            }
    }

    }
    
</style>