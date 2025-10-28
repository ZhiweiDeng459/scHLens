<template>
  <div class="container">
    <b class="title">Help</b>
    <p class="body">
      To help users learn how to use scHLens more quickly, we provide a usage
      guide here.
    </p>

    <p class="body">
      <b>Step 1.</b> Select a dataset.
      <router-link to="/tutorial/selectDataset"
        >[how to select dataset]</router-link
      >
    </p>

    <p class="body">
      <b>Step 2.</b> Config and run a pipeline.
      <router-link to="/tutorial/startaAnalysisPipeline"
        >[how to start a pipeline]</router-link
      >
    </p>

    <p class="body">
      <b>Step 3.</b> Explore your intersted pattern with the interactive views.
      <router-link to="/tutorial/analysiswithViews"
        >[how to analyze with views]</router-link
      >
    </p>

    <p class="body">
      <b>Step 4.</b> Annotate groups using manual labeling or our automated
      label type recommendation method.
      <a
        href="#"
        @click.prevent="
          tutorial_jump('/tutorial/analysiswithViews/2/8', 'tAnnotation')
        "
        title="how to annotate cells"
        >[how to annotate cells]</a
      >
    </p>

    <p class="body">
      <b>Step 5.</b>When you find that some cell groups may contain subgroups,
      you can select them using
      <a
        href="#"
        @click.prevent="
          tutorial_jump('/tutorial/analysiswithViews/2/1', 't_SelectCells')
        "
        title="how to select cells"
        >[how to select cells]</a
      >
      , configure a local pipeline, and then create a local plot. These
      subgroups will be displayed in Local Plot.
      <a
        href="#"
        @click.prevent="
          tutorial_jump(
            '/tutorial/hierarchicalExploration/1',
            't_CrateLocalPlot'
          )
        "
        title="how to create Local Plot"
        >[how to create Local Plot]</a
      >
    </p>

    <p class="body">
      <b>Step 6.</b> Once all annotations are complete, you can merge the
      annotation results and favorable visual patterns from the Local Plots into
      the Global Plot to derive the final annotation atlas.
      <a
        href="#"
        @click.prevent="
          tutorial_jump(
            '/tutorial/hierarchicalExploration/4',
            't_Merge'
          )
        "
        title="how to merge local results"
        >[how to merge local results</a
      >
    </p>

    <p class="body">
      You can learn scHLens by reproducing three examples:
      <router-link to="/tutorial/case1"
        >1. PBMC-3K</router-link
      >
      &nbsp; &nbsp;
      <router-link to="/tutorial/case2"
        >2. Human Pancreas</router-link
      >
      &nbsp; &nbsp;
      <router-link to="/tutorial/case3"
        >3. Human Lung Cancer</router-link
      >
    </p>

  </div>
</template>

<script>
import eventBus from "@/utils/eventBus.js";

export default {
  name: "THelp",
  methods: {
    tutorial_jump(index, id) {
      eventBus.$emit("TutorialJump", index, id);
    },
  },
};
</script>

<style lang="less" scoped>
.container {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  padding: 40px 100px;
  box-sizing: border-box;
  overflow: auto;
}

.title {
  font-size: 40px;
  margin: 40px 0px;
}

h1 {
  margin-top: 60px;
}

.h2 {
  margin-top: 40px;
  font-size: 29px;
}

.body {
  font-size: 22px;
  line-height: 200%;
  margin: 10px 0px;
}

.list {
  margin: 0px 20px;
}

.table-container {
  display: flex;
  flex-direction: column;
  align-items: center;
  font-size: 20px;

  table {
    border-collapse: collapse;
    /* 合并边框为单线 */
    border: 1px solid #333;
    /* 外框 */
  }

  th,
  td {
    border: 1px solid #333;
    /* 单元格框线 */
    padding: 6px 10px;
    /* 内边距 */
    text-align: center;
    /* 居中 */
  }
}

.image-container {
  display: flex;
  align-items: center;
  flex-direction: column;
  margin: 30px 0px;

  .image {
    border: 2px solid gray;
    box-shadow: 2px 2px 18px gray;
    max-width: 80%;
  }

  .image-title {
    font-size: 25px;
    margin: 15px 0px;
  }
}
</style>