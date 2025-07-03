<template>
  <div class="container">
    <b class="title">Start an Analysis Pipeline</b>
    <p class="body">
      After selecting the dataset, the next step is to configure and start the
      analysis pipeline. The pipeline is composed of seven modules that work
      together to process the data：
    </p>
    <ol>
      <li class="body">
        <b>Downsampling</b>: Downsample the dataset to speed up the pipeline and reduce
        memory usage. (Node: The "Stratified" method in the Downsampling module of the pipeline requires a dataset with
        labels. If you want to upload a dataset with labels, you should upload it in h5ad format and set the labels in
        AnnData.obs['label'].)
      </li>
      <li class="body">
        <b>Quality Control</b>: Filter low-quality cells and genes to ensure
        downstream analysis quality.
      </li>
      <li class="body">
        <b>Normalization</b>: Adjust the scRNA-seq counts to a common scale.
      </li>
      <li class="body">
        <b>Gene Selection</b>: Filter uninformative or noise genes and only
        Focus on the biological variable genes.
      </li>
      <li class="body">
        <b>Visualization</b>: Embed high-dimensional cell expression data into
        two-dimensional space, allowing scRNA-seq gene expression information to
        be visualized at the single-cell granularity.
      </li>
      <li class="body">
        <b>Clustering</b>: Group cells with similar expression profiles, further
        assisting analysts in identifying cell types.
      </li>
      <li class="body">
        <b>Differential Expressed Gene Identification (DEG Identification)</b>: Discover differential expressed genes
        that characterize
        each cell group.
      </li>
    </ol>
    <div class="image-container">
      <img class="image" src="tutorials/img/Pipeline.jpg" />
      <b class="image-title">Figure 1. Pipeline Interface</b>
    </div>
    <p class="body">
      To configure the pipeline, click the blue "Set Pipeline" button in the
      top-left corner (The red highlighted area in Figure 1) to open the
      pipeline dialog.
    </p>
    <p class="body">
      Within the pipeline dialog, individual module can be added by clicking the
      "Add" button, all modules can be added at once with the "All" button, and
      the added modules can be cleared with the "Clear" button (Area A in Figure
      1).
    </p>
    <p class="body">
      The dialog will display the added modules and allow users to configure the
      settings for each module (Area B in Figure 1).
      <b>In these settings, if you do not wish to configure a specific
        parameter, you should leave its input field blank.</b>
    </p>
    <p class="body">
      After completing the pipeline configuration, click the green "Run" button
      (Area C in Figure 1) to start the pipeline.
    </p>
  </div>
</template>

<script>
export default {
  name: "TStartaAnalysisPipeline",
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
  font-size: 29px
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