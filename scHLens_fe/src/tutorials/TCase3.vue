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
    <h1 id="t_re_case3">
      Case3: scHLens hierarchically identifies detailed cell subtypes from the large-scale human lung cancer dataset.
    </h1>
    <p class="body">
      The target dataset for this case is Human Lung Cancer dataset, containing 256,379 cells with 32,805 genes. The
      original annotation of this dataset contains 9 types, including Airway Epithelium, B cells, CD45+, Endothelial,
      Epithelial, Fibroblasts, Granulocytes, Myeloid, and T cells.
    </p>
    <p class="body">
      We list the essential pipeline parameters for each plot as follows. We list the essential pipeline parameters for
      each plot as follows. Users can reproduce the results step by step following the descriptions in the manuscript.
      <!--Global Plot 1 Parameter-->
    <div class="table-container">
      <table>
        <tr>
          <th>Module</th>
          <th>Method</th>
          <th>Parameter</th>
          <th>Value</th>
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
          <td>3000</td>
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
          <td>0.6</td>
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
          <td>15</td>
        </tr>
        <tr>
          <td>Clustering</td>
          <td>leiden</td>
          <td>resolution</td>
          <td>0.015</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>wilcoxon-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-1: The pipeline parameter of Global Plot &#9312; in Case 3
      </figcaption>
    </div>



    <!--Local Plot 2 Parameter-->
    <div class="table-container">
      <table>
        <tr>
          <th>Module</th>
          <th>Method</th>
          <th>Parameter</th>
          <th>Value</th>
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
          <td>0.1</td>
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
          <td>0.05</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>wilcoxon-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-2: The pipeline parameter of Local Plot &#9313; in Case 3
      </figcaption>
    </div>



    <!--Local Plot 3 Parameter-->
    <div class="table-container">
      <table>
        <tr>
          <th>Module</th>
          <th>Method</th>
          <th>Parameter</th>
          <th>Value</th>
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
          <td>Gene Selection</td>
          <td>HighlyVariable</td>
          <td>top Genes</td>
          <td>1000</td>
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
          <td>0.1</td>
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
          <td>0.15</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>wilcoxon-test</td>
          <td>nGenes</td>
          <td>100</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-3: The pipeline parameter of Local Plot &#9314; in Case 3
      </figcaption>
    </div>

    <!--Local Plot 4 Parameter-->
    <div class="table-container">
      <table>
        <tr>
          <th>Module</th>
          <th>Method</th>
          <th>Parameter</th>
          <th>Value</th>
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
          <td>scry</td>
          <td>top Genes</td>
          <td>300</td>
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
          <td>0.6</td>
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
          <td>0.05</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>t-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-4: The pipeline parameter of Local Plot &#9315; in Case 3
      </figcaption>
    </div>

    <!--Local Plot 5 Parameter-->
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
          <td></td>
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
          <td>Yes</td>
        </tr>
        <tr>
          <td>Gene Selection</td>
          <td>HighlyVariable</td>
          <td>top Genes</td>
          <td>2000</td>
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
          <td>15</td>
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
          <td>0.1</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>t-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-5: The pipeline parameter of Local Plot &#9316; in Case 3
      </figcaption>
    </div>


    <!--Local Plot 6 Parameter-->
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
          <td></td>
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
          <td>Yes</td>
        </tr>
        <tr>
          <td>Gene Selection</td>
          <td>scry</td>
          <td>top Genes</td>
          <td>1000</td>
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
          <td>15</td>
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
          <td>0.1</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>t-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-6: The pipeline parameter of Local Plot &#9317; in Case 3
      </figcaption>
    </div>


    <!--Local Plot 7 Parameter-->
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
          <td></td>
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
          <td>Yes</td>
        </tr>
        <tr>
          <td>Gene Selection</td>
          <td>scry</td>
          <td>top Genes</td>
          <td>500</td>
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
          <td>0.1</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>t-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-7: The pipeline parameter of Local Plot &#9318; in Case 3
      </figcaption>
    </div>

    <!--Local Plot 8 Parameter-->

    <div class="table-container">
      <table>
        <tr>
          <th>Module</th>
          <th>Method</th>
          <th>Parameter</th>
          <th>Value</th>
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
          <td>scry</td>
          <td>top Genes</td>
          <td>300</td>
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
          <td>0.6</td>
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
          <td>0.05</td>
        </tr>
        <tr>
          <td>DEG Identification</td>
          <td>wilcoxon-test</td>
          <td>nGenes</td>
          <td>50</td>
        </tr>
        <tr>
          <td>Cell Communication</td>
          <td>CellChat</td>
          <td>Organism</td>
          <td>Human</td>
        </tr>
      </table>
      <figcaption style="margin-top: 6px; font-weight: bold">
        Table 3-8: The pipeline parameter of Local Plot &#9319; in Case 3
      </figcaption>
    </div>


  </div>
</template>

<script>
export default {
  name: "TCase3",
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
  margin-bottom: 30px;

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