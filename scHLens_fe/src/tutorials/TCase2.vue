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
    <h1 id="t_re_case2">
      Case2: scHLens identifies endocrine subtypes and rare cell types in the
      human pancreatic dataset
    </h1>
    <p class="body">
      The target dataset for this case is Human Pancreatic dataset, which is
      representative dataset of human pancreatic cell atlas, obtained from an
      automated platform that utilizes FACS and CEL-Seq2. It was originally
      annotated as nine major cell types, including five kinds of Endocrine
      cells (Alpha cells, Beta cells, Delta cells, Epsilon cells and Gamma
      cells), two kinds of Exocrine cells (Ductal cells and Acinar cells),
      Mesenchymal cells and Endothelial cells.
    </p>
    <p class="body">
      <b>Step 1.</b>Select the Human Pancreatic dataset (sample) and then Set the global pipeline as Table 2-1 and then
      obtain the unannotated Global Plot &#9312; (Shown
      as Figure 2-1).
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
          <td>T-SNE</td>
          <td>pre-DR</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>T-SNE</td>
          <td>perplexity</td>
          <td>30</td>
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
          <td>0.03</td>
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
        Table 2-1: The pipeline parameter of Global Plot &#9312; in Case 2
      </figcaption>
    </div>
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure2-1.jpg" />
      <b class="image-title">Figure 2-1</b>
    </div>
    <p class="body">
      <b>Step 2.</b>The unannotated Global Plot &#9312; consists six
      groups (Shown as Figure 2-1). We selected Human - PanglaoDB Augmented 2021 as the gene set
      database and Enrichr as the annotation recommendation method to annotate
      each group. During the annotation process, we applied different
      logfoldchange thresholds to filter marker genes for each group: c_0: 3,
      c_1: 2, c_2: 5, c_3: 4, c_4: 2, c_5: 2. According to the recommendations
      from the Annotation module, c_2 - c_5 were accurately
      identified: Acinar cells, Ductal cells, Mesenchymal cells, and Endothelial
      cells. The clusters c_0 and c_1 exhibit relatively clear several blocks,
      and therefore, they were selected for more detailed clustering (Figure 2-2).
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure2-2.jpg" />
      <b class="image-title">Figure 2-2</b>
    </div>
    <p class="body">
      <b>Step 3.</b>First, we selected c_0 and created a local pipeline as shown in parameter
      Table 2-2, resulting in Local Plot &#9313;. In Local Plot &#9313;, there are three
      groups, and we again used Human - PanglaoDB Augmented 2021 as the gene set
      database and Enrichr as the annotation recommendation method to annotate
      each group. We applied different logfoldchange thresholds to filter marker genes for each group: c_0: 1,
      c_1: 1, c_2: 3. The cell type recommendations are shown as Figure 2-3 (B). The final cell type was selected as the
      most significant recommended type or the parent type of the top recommended types. Therefore, we can annotate it as
      Figure 2-3 (A).
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
          <td>T-SNE</td>
          <td>pre-DR</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>T-SNE</td>
          <td>perplexity</td>
          <td>30</td>
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
          <td>0.03</td>
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
        Table 2-2:  The pipeline parameter of Local Plot &#9313; in Case 2
      </figcaption>
    </div>
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure2-3.jpg" />
      <b class="image-title">Figure 2-3</b>
    </div>
    <p class="body">
      <b>Step 3.</b> To further validate our annotation results, we examined the expression patterns of marker genes
      corresponding to the annotated cell types. By entering the marker genes into the gene search box and adding them to
      the "Current Gene" list, we visualized their expression in the Gene Projection View. As shown in Figure 2-3 (C), the
      marker genes of Alpha Cells (GCG, TTR, MAFB, FAP, GLS) were highly expressed in cluster c_0, confirming the
      correctness of the annotation. Similarly, the annotations of c_1 and c_2 as Gamma Cells and Myeloid Cells,
      respectively, were also accurate. The annotation in the the original publication of c_2 was Ductal Cell; however,
      the marker genes of
      Ductal Cells were not highly expressed in c_2, indicating a misannotation.
    </p>

    <p class="body">
      <b>Step 4.</b>Returning to Global Plot &#9312;, we selected c_1 and created a local pipeline as shown in Table 2-3, resulting
      in Local Plot &#9314;.In Local Plot &#9314;, there are two
      groups, and we again used Human - PanglaoDB Augmented 2021 as the gene set
      database and Enrichr as the annotation recommendation method to annotate
      each group. We applied same logfoldchange thresholds to filter marker genes for each group: c_0: 2,
      c_1: 2. The final cell type was selected as the
      most significant recommended type or the parent type of the top recommended types. Therefore, we can annotate it as
      Figure 2-4.
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
          <td>500</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>T-SNE</td>
          <td>pre-DR</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>T-SNE</td>
          <td>perplexity</td>
          <td>30</td>
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
        Table 2-3 The pipeline parameter of Local Plot &#9315; in Case 2
      </figcaption>

      <div class="image-container">
        <img class="image" src="tutorials/img/repo-Figure2-4.jpg" />
        <b class="image-title">Figure 2-4</b>
      </div>
    </div>

    <p class="body">
      We list the essential pipeline parameters for remaining plots:
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
          <td>80</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>T-SNE</td>
          <td>pre-DR</td>
          <td>No</td>
        </tr>
        <tr>
          <td>Visualization</td>
          <td>T-SNE</td>
          <td>perplexity</td>
          <td>30</td>
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
        Table 2-4: The parameter of pipeline in Case &#9316; in Case 2
      </figcaption>
    </div>


  </div>
</template>

<script>
export default {
  name: "TCase2",
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