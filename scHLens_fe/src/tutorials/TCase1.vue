<template>
  <div class="container">
    <b class="title">Case reproduction</b>
    <p class="body">
      In our manuscript, we conducted three case studies to demonstrate the
      effectiveness of scHLens in real-world analysis scenarios. These cases are
      as follows: Case 1, clustering and annotating the PBMC-3K dataset; Case 2,
      scHLens identifies endocrine subtypes and rare cell types in the human
      pancreatic dataset; and Case 3, scHLens hierarchically identifies detailed
      cell subtypes from the large-scale human lung cancer dataset.
    </p>
    <p class="body">
      We will next describe the analysis process for each case in detail,
      allowing users to reproduce our results in manuscript. <b>However, some
        methods in our pipeline involve stochastic processes, so the generated
        visualizations may not be exactly identical to those shown in the
        manuscript figures (for example, the relative positions of groups in the
        Cell Projection View). We have tested and confirmed that such randomness
        only introduces minor variations and does not affect the overall
        conclusions.</b>
    </p>
    <h1 id="t_re_case1">
      Case1: clustering and annotating the PBMC-3K dataset
    </h1>
    <p class="body">
      The target dataset for this case is PBMC-3K, which consists of peripheral
      blood mononuclear cells (PBMCs) from a healthy donor. It represents a
      diverse set of immune cell types, including T cells, B cells, natural
      killer (NK) cells, and monocytes.
    </p>
    <p class="body">
      <b>Step 1.</b>
      Select the PBMC-3K (sample) dataset and then Set the global pipeline as Table 1 and then obtain the unannotated plot
      (Shown
      as Figure 1-1).
    </p>
    <div class="table-container">
      <table>
        <tr>
          <th>Module</th>
          <th>Method</th>
          <th>Parameter</th>
          <th>Value</th>
        </tr>
        <tr>
          <td>Quality Control</td>
          <td>Filter Outlying Cells</td>
          <td>min Genes</td>
          <td>200</td>
        </tr>
        <tr>
          <td>Quality Control</td>
          <td>Filter Outlying Cells</td>
          <td>max Genes</td>
          <td>2500</td>
        </tr>
        <tr>
          <td>Quality Control</td>
          <td>Filter Outlying Genes</td>
          <td>min Cells</td>
          <td>3</td>
        </tr>
        <tr>
          <td>Quality Control</td>
          <td>Filter Outlying Genes</td>
          <td>max Cells</td>
          <td></td>
        </tr>
        <tr>
          <td>Quality Control</td>
          <td>MT-Gene Control</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
        <tr>
          <td>Quality Control</td>
          <td>MT-Gene Control</td>
          <td>pct Counts</td>
          <td>5</td>
        </tr>
        <tr>
          <td>Normalization</td>
          <td>CPM Normalize</td>
          <td>-</td>
          <td>Yes</td>
        </tr>
        <tr>
          <td>Normalization</td>
          <td>Logarithmize</td>
          <td>-</td>
          <td>Yes</td>
        </tr>
        <tr>
          <td>Normalization</td>
          <td>Regress</td>
          <td>-</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Normalization</td>
          <td>Scale</td>
          <td>-</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Gene Selection</td>
          <td>HighlyVariable</td>
          <td>top Genes</td>
          <td>1500</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>UMAP</td>
          <td>pre-DR</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>UMAP</td>
          <td>min Dist</td>
          <td>0.5</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>UMAP</td>
          <td>nNeighbors</td>
          <td>10</td>
        </tr>
        <tr>
          <td>Clustering</td>
          <td>leiden</td>
          <td>nNeighbors</td>
          <td>10</td>
        </tr>
        <tr>
          <td>Clustering</td>
          <td>leiden</td>
          <td>resolution</td>
          <td>0.5</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>t-test</td>
          <td>nGenes</td>
          <td>100</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 1
      </figcaption>
    </div>
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure1-1.jpg" />
      <b class="image-title">Figure 1-1</b>
    </div>
    <p class="body">
      <b>Step 2.</b>
      After executing the above pipeline, we obtained a plot consisting of eight
      groups (Figure 1-1). We then opened the Group Panel, selected Human -
      Azimuth 2023 as the gene set database and Enrichr as the annotation
      recommendation method to annotate each group. During the annotation
      process, we applied different logfoldchange thresholds to filter marker
      genes for each group: c_0: 1, c_1: 4, c_2: 4, c_3: 2.3, c_4: 3, c_5: 4,
      c_6: 3, and c_7: 5. Thie cell type recommendation results are shown in Figure 1-2.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure1-2.jpg" />
      <b class="image-title">Figure 1-2</b>
    </div>
    <p class="body">
      <b>Step 3.</b>
      The final cell type was selected as the most significant recommended type or the parent type of the top
      recommended types. The final annotation, shown in Figure 1-3 of the manuscript, is similar to the
      original annotation presented in Figure 1-1.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure1-3.jpg" />
      <b class="image-title">Figure 1-3</b>
    </div>

  </div>
</template>

<script>
export default {
  name: "TCase1",
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
}</style>