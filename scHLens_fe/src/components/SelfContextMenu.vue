<template>
    <div class="menu-container"
        ref="menu-container"
        tabindex="0"
        v-show="visibility">
        <div v-for="item in items"
            class="menu-item"
            :key="item.name"
            @click="item.callback"
            >
            <img class="item-icon" :src="item.icon"/>
            <a class="item-text">{{item.name}}</a>
        </div>
    </div>

</template>

<script>

import * as d3 from "d3";


export default {
    name:"SelfContextMenu",
    props:[
        'items', // {'name':...,callback:...,icon:...}
        '_mounted', //(this,menuDOM,base)=>{}
    ],
    data(){
        return {
          visibility:false  
        }
    },
    methods:{
        setVisibility(isVis){
            this.visibility = isVis
        },
        getVisibility(){
            return this.visibility
        },

    },
    mounted(){
        
        const menuDOM = this.$refs['menu-container']
        //找到父节点
        const base = menuDOM.parentNode
        //找到根节点
        const root = document.getElementById("app")

        //_外部函数
        if(typeof this._mounted == 'function')
            this._mounted(this,menuDOM,base);//当前vue对象，容器DOM，父节点DOM

        document.addEventListener('click',(e)=>{ //任意的单击事件都会关闭菜单栏
            this.setVisibility(false)
        })

        document.addEventListener('contextmenu',(e)=>{

            if(e.target === menuDOM || menuDOM.contains(e.target)){
                e.preventDefault();
                this.setVisibility(false)
            }
            else if(base.contains(e.target)){
                e.preventDefault();
                //显示菜单
                this.setVisibility(true)
                //调整位置
                menuDOM.style.top = e.layerY + 'px'
                menuDOM.style.left = e.layerX + 'px' 
            }
            else{
                e.preventDefault();
                this.setVisibility(false)
            }
        })
    



        d3.select(this.$refs['menu-container'])
            .selectAll('.menu-item') //菜单项相关
            .data(this.items)
            .on('mouseenter.contextmenu',function(e,d){
                this.classList.add('selected')
                this.classList.remove('unselected')
            })
            .on('mouseleave.contextmenu',function(e,d){
                this.classList.add('unselected')
                this.classList.remove('selected')
            })
    }
}


</script>

<style lang="less">
    .menu-container{
        width: 280px;
        border:1px solid lightgray;
        border-radius: 2px;
        position: absolute;
        z-index: 1000;
        box-shadow: 0px 0px 2px lightgray;
        background-color:#eeeeee;
        .menu-item{
            height: 30px; //菜单项高度
            margin: 4px 0px;
            // border:1px solid lightgray;
            // background-color: white;
            cursor: pointer;
            display: flex;
            align-items: center;

            .item-icon{
                flex:0 0 45px;
                margin:0px 2px 0px 5px;
                height: 25px;
            }
            .item-text{
                font-size:18px;
                font-family:YaHei;
            }
        }

        .unselected{
            
        }
        .selected{
            background-color: lightblue;
        }

    }
</style>