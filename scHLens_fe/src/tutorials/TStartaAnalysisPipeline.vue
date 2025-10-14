<template>
  <div class="container">
    <b class="title">Start an Analysis Pipeline</b>

    <h1 id="t_start_a_pipeline">How to start a pipeline</h1>
    <p class="body">
      After selecting the dataset, the next step is to configure and start the
      analysis pipeline. The pipeline is composed of seven modules that work
      together to process the data：
    </p>
    <ol>
      <li class="body">
        <b>Downsampling</b>: Downsample the dataset to speed up the pipeline and
        reduce memory usage. (Node: The "Stratified" method in the Downsampling
        module of the pipeline requires a dataset with labels. If you want to
        upload a dataset with labels, you should upload it in h5ad format and
        set the labels in AnnData.obs['label'].)
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
        <b>Differential Expressed Gene Identification (DEG Identification)</b>:
        Discover differential expressed genes that characterize each cell group.
      </li>
      <li class="body">
        <b>Cell Communication</b>:
        Detect which cell types send and receive signals (ligand - receptor pairs) and Help to map cellular interaction networks.
      </li>
    </ol>
    <div class="image-container">
      <img class="image" src="tutorials/img/Pipeline.jpg" />
      <b class="image-title">Figure 1. Pipeline Interface</b>
    </div>
    <p class="body">
      To configure the pipeline, click the blue "Set Pipeline" button in the
      left side of the system (The red highlighted area in Figure 1) to open the
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
      <b
        >In these settings, if you do not wish to configure a specific
        parameter, you should leave its input field blank.</b
      >
    </p>
    <p class="body">
      After completing the pipeline configuration, click the green "Run" button
      (Area C in Figure 1) to start the pipeline.
    </p>
    <h1 id="t_set_parameters">How to set parameters</h1>
    <p class="body">
      We provided detailed explanations for the less intuitive parameters in the
      module’s methods, helping users understand their functions and effects,
      and recommended suitable value ranges.
    </p>
    <p id="t_tSNE_perplexity" class="body">
      <b>Visualization - t-SNE - perplexity:</b> In t-SNE, perplexity can be
      regarded as a smooth measure of the effective number of nearest neighbors
      considered when computing pairwise similarities in the high-dimensional
      space. On the one hand, with a fixed dataset size, a smaller perplexity
      emphasizes local structure, causing the algorithm to capture fine-grained
      relationships and often yielding more fragmented visual clusters in the
      embedding. In contrast, a larger perplexity emphasizes global structure,
      promoting the preservation of broader neighborhood relationships and
      producing fewer, more cohesive visual clusters, albeit sometimes at the
      expense of local detail. On the other hand, as dataset size increases,
      each point has more valid neighbors, necessitating a larger perplexity.
      For small datasets (fewer than 1,000 cells), a perplexity of 10–40 is
      recommended; for larger datasets (more than 1,000 cells), a perplexity of
      20–50 is more appropriate. In scHLens, we set the default perplexity to
      30, which is suitable for most scenarios.
    </p>
    <p id="t_UMAP_minDist" class="body">
      <b>Visualization - UMAP - min Dist:</b> In UMAP, the min Dist parameter
      controls the minimum distance between points in the low-dimensional
      embedding. It essentially regulates how tightly points are packed together
      in the resulting visualization. A smaller min Dist value allows points to
      cluster more closely, highlighting local structure and fine-grained
      patterns, which makes clusters appear visually denser. Conversely, a
      larger min Dist spreads points farther apart, emphasizing global structure
      and making clusters more visually separated, often resulting in a more
      uniform and less crowded layout. For visualizing scRNA-seq subpopulations,
      a recommended min Dist range is 0–0.5.
    </p>
    <p id="t_UMAP_nNeighbors" class="body">
      <b>Visualization - UMAP - nNeighbors:</b> In UMAP, the nNeighbors
      parameter determines the number of nearest neighbors considered when
      constructing the local neighborhood graph. It directly influences the
      balance between capturing local versus global structure in the embedding.
      On the one hand, with a fixed dataset size, a smaller nNeighbors value
      focuses more on local relationships, which can reveal fine-grained visual
      clusters and local patterns but may fragment global structure. Conversely,
      a larger nNeighbors considers more points for each neighborhood,
      emphasizing the global structure and producing smoother, more connected
      visual clusters, though at the cost of local detail. On the other hand, as
      dataset size increases, each point has more valid neighbors, necessitating
      a larger nNeighbors. For small datasets (fewer than 1,000 cells), a
      nNeighbors of 5–20 is recommended; for larger datasets (more than 1,000
      cells), a nNeighbors of 10–50 is more appropriate. In scHLens, we set the
      default nNeighbors to 15, which is suitable for most scenarios.
    </p>
    <p id="t_Leiden_nNeighbors" class="body">
      <b>Clustering - Leiden - nNeighbors:</b> In the Leiden algorithm, the
      nNeighbors parameter determines the number of nearest neighbors used when
      constructing the k-nearest neighbor (kNN) graph, which serves as the
      foundation for community detection. On the one hand, with a fixed dataset
      size, a smaller nNeighbors value restricts the neighborhood, making the
      algorithm more sensitive to fine-grained local structures and often
      leading to the identification of smaller, more fragmented clusters.
      Conversely, a larger nNeighbors value expands the neighborhood,
      emphasizing global connectivity and resulting in fewer, larger clusters,
      though potentially at the expense of local detail. On the other hand, as
      dataset size increases, each point has more valid neighbors, necessitating
      a larger nNeighbors. For small datasets (fewer than 1,000 cells), a
      nNeighbors of 5–20 is recommended; for larger datasets (more than 1,000
      cells), a nNeighbors of 10–50 is more appropriate. In scHLens, we set the
      default nNeighbors to 10, which is suitable for most scenarios.
    </p>
    <p id="t_Leiden_resolution" class="body">
      <b>Clustering - Leiden - resolution:</b>In the Leiden algorithm, the
      resolution parameter controls the granularity of the detected clusters. A
      higher resolution value leads to a larger number of smaller, more
      fine-grained clusters, effectively splitting communities into finer
      subgroups. Conversely, a lower resolution produces fewer, larger clusters,
      capturing broader groupings at the cost of potentially merging distinct
      subpopulations. We recommend setting the resolution parameter within the
      range of 0 to 2.
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