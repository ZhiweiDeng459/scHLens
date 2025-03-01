<template>
  
    <div class="tutorial-container">

        <router-view/>

    </div>

</template>

<script>
export default {
    name: "TutorialContainer",

}
</script>

<style>
    .tutorial-container{
        width:100%;
        height: 100%;
        /* overflow: auto; */
    }
</style>