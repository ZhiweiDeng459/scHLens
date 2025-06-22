<template>
    <div>
        <div class="state-panel-container">
            <div id="state">
                <img id="state-icon" src="icons\finish.svg"/>
                <a id="state-text" style="color:white">loading....</a>
            </div>
        </div>
    </div>
</template>

<script>
export default {
    name:"StatePanel",
    data(){
        return {

        }
    },
    methods:{

    },
}
</script>

<style lang="less">
.state-panel-container{
    width:100%;
    height:100%;
    display: flex;
    #state{
        width:100%;
        height:100%;
        display: flex;
        align-items: center;
        #state-icon{
            height:20px;
            margin:0px 10px 0px 0px;
        }
        #state-text{
            font-size:20px
        }
    }
}
</style>