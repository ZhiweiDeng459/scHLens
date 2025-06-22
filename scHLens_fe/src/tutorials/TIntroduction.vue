<template>
  <div class="container">
    <b class="title">Introduction</b>
    <p class="body">
      scHLens is an interactive analysis platform for hierarchically and
      interactively exploring cellular heterogeneity of scRNA-seq dataset, with
      various visualization views. scHLens focused on the cell type annotation
      of scRNA-seq. The main novelty of scHLens includes two aspects. Firstly,
      for each layer, scHLens provided a customized pipeline for gene selection,
      cell visualization, clustering, differential expressed gene analysis and
      cell type annotation. Then, scHLens provided a hierarchical analytical
      mode to explore cellular heterogeneity, and proposed a novel
      hierarchically visualization method to integrate different levels of cells
      subgroups. This web platform allows users to easily select, export, delete
      and reanalyze specific cell groups of their interest.
    </p>
    <p class="body">
      The system's workflow is composed of five steps: <b>Enter a Job</b>,
      <b>Select a Dataset</b>, <b>Start the User-defined Pipeline</b>,
      <b>Analyze with Interactive Views</b>, <b>Hierarchically Explore</b>. Each
      step will be explained in detail below.
    </p>
    <p class="body">
      <b
        >To better demonstrate the cases in our manuscript, we have integrated
        predefined demos into scHLens. For details on how to use them, please
        refer to "Enter a Job - Create a New Job & Load an Existing Job & Use
        Samples - Use Samples".</b
      >
    </p>
  </div>
</template>

<script>
export default {
  name: "TIntroduction",
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

h2 {
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