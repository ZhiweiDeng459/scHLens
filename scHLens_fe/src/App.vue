<template>
    <div id="app">
        <div style="display:flex">
            <div class="TitleBoard">
                <b class="SystemTitle">scHLens</b>
            </div>
            <el-menu
                style="flex:1 1 0"
                :default-active="$route.path"  
                mode="horizontal"
                background-color="#24292f"
                text-color="#fff"
                active-text-color="#90e36b"
                router>
                <el-menu-item index="/home" >
                    <a class="menu-item-text" slot="title">Home</a>
                </el-menu-item>
                <el-menu-item index="/system">
                    <span class="menu-item-text" slot="title">System</span>
                </el-menu-item>
                <el-menu-item index="/tutorial">
                    <span class="menu-item-text" slot="title">Tutorial</span>
                </el-menu-item>
                <el-menu-item index="/contact">
                    <span class="menu-item-text" slot="title">Contact us</span>
                </el-menu-item>
            </el-menu>
        </div>
        <keep-alive>
            <router-view/>
        </keep-alive>

        <!--游走小信息板-->
        <InfoPanel ref="InfoPanel"/>

        <!--固定大信息板-->
        <MessageBoard ref="MessageBoard"/>
        

    </div>
</template>

<script>
import Vue from "vue";
import {Step,Steps,Message,Menu, MenuItem, Submenu,Notification} from "element-ui";
import InfoPanel from "./components/InfoPanel.vue";
import MessageBoard from "./components/MessageBoard.vue";
import AnnoRecomAddiInfoPanel from "@/components/AnnoRecomAddiInfoPanel";
import {InstanceClose} from '@/utils/interface';
import axios from "axios";
import Tutorial from "@/Tutorial";
Vue.component(Menu.name,Menu)
Vue.component(MenuItem.name,MenuItem)
Vue.component(Submenu.name,Submenu)
Vue.prototype.$message = Message;
Vue.prototype.$notify = Notification;

export default {
    name: "App",
    components: {
        InfoPanel,
        MessageBoard,
    },
    computed: {
        curData() {
            return this.$store.state.curData;
        },
        dataList(){
            return this.$store.state.dataList;
        },
        pipelineEntityArray(){
            return this.$store.state.pipelineEntityArray;
        },
        JobId(){
            return this.$store.state.JobId
        },
        socketIns(){
            return this.$store.state.socketIns;
        }

    },
    methods:{
        closeHandler(event){
            //网页关闭或者刷新
            InstanceClose()
        },
    },

    data(){
        return{

        }
    },


    created(){
        //绑定网页关闭处理事件
        window.addEventListener('beforeunload',this.closeHandler)
        //设置axios的baseURL
    },

    mounted(){
        //set InfoPanel
        this.$store.commit("setInfoPanel", this.$refs['InfoPanel']);
        //set MessageBoard
        this.$store.commit("setMessageBoard", this.$refs['MessageBoard']);
        // //set AnnoRecomAddiInfoPanel
        // this.$store.commit("setAnnoRecomAddiInfoPanel", this.$refs['AnnoRecomAddiInfoPanel']);


        //预渲染Tutorial
        const initialRoute = this.$router.currentRoute.fullPath;
            // 临时切换到 /tutorial
            this.$router.push("/tutorial").then(() => {
            // 预加载完成后再切回初始路由
            this.$router.push(initialRoute);
        });

        
    }

};
</script>

<style scoped lang="less">
#app {
    font-family: Avenir, Helvetica, Arial, sans-serif;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
    height: 100vh;
    width: 100vw;
    overflow: hidden;
    display: flex;
    flex-direction: column;
    position: relative;
    .TitleBoard{

        background-color:#24292f;
        height: 60px;
        width:  200px;
        display: flex;
        align-items: center;
        .SystemTitle{
            font-size: 40px;
            margin-left : 35px;
            font-family: AslinaBold;
            color: #90e36b;
            margin-top:15px;
            margin-bottom:20px;
        }
    }
    
    .menu-item-text{
        font-size: 20px;
        padding: 5px;
        
    }

}


</style>
