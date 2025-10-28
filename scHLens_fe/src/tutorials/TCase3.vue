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
      The target dataset for this case is Human Lung Cancer dataset. which is from 9 treatment-naïve LUAD patients from
      Brigham and Women’s Hospital. This dataset can be accessed at
      <a herf="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013"
        target="_blank">https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013.</a>
      This dataset contains 256,379 cells with 32,805 genes.
      The original annotation of this dataset contains 9 types, including Airway Epithelium, B cells, CD45+, Endothelial,
      Epithelial, Fibroblasts, Granulocytes, Myeloid, and T cells.
    </p>


    <!-- step 1-->
    <p class="body">
      <b>Step 1.</b>Select the Human Pancreatic dataset (sample) and then Set
      the global pipeline as Table 3-1 and then obtain the unannotated Global
      Plot &#9312; (Shown as Figure 3-1).
    </p>

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

    <!--Figure 3-1-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-1.png" />
      <b class="image-title">Figure 3-1</b>
    </div>



    <p class="body">
      <b>Step 2.</b>The unannotated Global Plot &#9312; consists ten groups
      (Shown as Figure 3-1). We selected Human - PanglaoDB Augmented 2021 as the
      gene set database and Enrichr as the annotation recommendation method to
      annotate each group. During the annotation process, we applied the same
      logfoldchange thresholds to filter marker genes for each group: logfoldchanges > 2.
      Based on the recommendations from the Annotation module, four cell types were accurately identified: B cells, Mast
      cells, Plasmacytoid Dendritic cells, Endothelial cells, and a subset of Epithelial cells. Clusters c_0, c_1, c_2,
      and c_3 exhibit distinct subtypes. (Figure 3-2) Therefore, we selected these four clusters along with all epithelial
      cell types
      for reanalysis.
    </p>

    <!--Figure 3-2-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-2.jpg" />
      <b class="image-title">Figure 3-2</b>
    </div>

    <p class="body">
      <b>Step 3.</b>First, we select cluster c_0 in Global Plot &#9312; (Figure 3-2) to run local pipeline as shown in
      Table 3-2.
      resulting in Local Plot &#9313;. In Local Plot
      &#9313;, there are three groups, and we again used Human - PanglaoDB
      Augmented 2021 as the gene set database and Enrichr as the annotation
      recommendation method to annotate each group. We applied different
      logfoldchange thresholds to filter marker genes for each group: c_0: 1,
      c_1: 3, c_2: 2. The cell type recommendations are shown as Figure 3-3.
      The final cell type was selected as the most significant recommended type
      or the parent type of the top recommended types. Therefore, we can
      annotate it as Figure 3-4. As shown in Figure 3-5, we utilized 2d-density plots to exhibit the marker genes for T
      cells and NK cells The higher expression of <em>CD3</em> in T cells and <em>FCGR3A</em> in NK cells, along with the lower
      expression of <em>CD3</em> in NK cells, supported these annotations.
    </p>

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

    <!--Figure 3-3-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-3.jpg" />
      <b class="image-title">Figure 3-3</b>
    </div>

    <!--Figure 3-4-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-4.jpg" />
      <b class="image-title">Figure 3-4</b>
    </div>

    <!--Figure 3-5-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-5.jpg" />
      <b class="image-title">Figure 3-5</b>
    </div>

    <p class="body">
      <b>Step 4.</b>We then performed a further analysis for sub-cluster c_2 in Local Plot &#9313; to create new Local
      Plot, &#9314; (The pipeline parameter is in Table 3-3). This is because we believed it contained multiple
      subtypes according to its EnrichR result (Figure 3-3). We identified three sub-clusters. We used Human - PanglaoDB
      Augmented 2021 as the gene set database and Enrichr as the annotation
      recommendation method to annotate each group. We applied different
      logfoldchange thresholds to filter marker genes for each group: c_0: 2,
      c_1: 2, c_2: 4. The cell type recommendations are shown as Figure 3-6.
      In Local Plot &#9314;. These sub-clusers were annotated as T
      cells, B cells, and B cells (Figure 3-7). The expression of markers (<em>CD3</em> for T cells;
      <em>CD19</em>, <em>CD79A</em>, <em>CD79B</em> for B cells) in Figure 3-8 validated the reliability of our annotations
    </p>

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

    <!--Figure 3-6-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-6.jpg" />
      <b class="image-title">Figure 3-6</b>
    </div>

    <!--Figure 3-7-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-7.jpg" />
      <b class="image-title">Figure 3-7</b>
    </div>

    <!--Figure 3-8-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-8.jpg" />
      <b class="image-title">Figure 3-8</b>
    </div>

    <p class="body">
      <b>Step 5.</b> Returning to the Global Plot &#9312;, we next select c_1 to generate a new Local Plot &#9315; (the
      pipeline
      parameters are listed in Table 3–4).
      In Local Plot &#9315;, we identified two sub-clusters. We used Human - PanglaoDB
      Augmented 2021 as the gene set database and Enrichr as the annotation
      recommendation method to annotate each group. We applied
      logfoldchange thresholds to filter marker genes for each group: c_0: 2,
      c_1: 2. The cell type recommendations are shown as Figure 3-9. Therefore, these sub-clusers were annotated as
      Dendritic Cells and Monocytes (Figure 3-10).
    </p>
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

    <!--Figure 3-9-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-9.jpg" />
      <b class="image-title">Figure 3-9</b>
    </div>

    <!--Figure 3-10-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-10.jpg" />
      <b class="image-title">Figure 3-10</b>
    </div>


    <p class="body">
      <b>Step 6.</b> Returning to the Global Plot &#9312;, according to the high expression of <em>EPCAM</em>, we select all
      epithelial cells (Figure 3-11) to create Local Plot &#9316; (the pipeline parameters are listed in Table 3–5).
      The Global Plot &#9316; was identified as nine sub-clussters. We used Human - PanglaoDB
      Augmented 2021 as the gene set database and Enrichr as the annotation
      recommendation method. We applied logfoldchange thresholds to filter marker genes: c_5: 4,
      c_7: 4, c_8: 4. The cell type recommendations are shown as Figure 3-12. Therefore, c_5, c_7, c_8 were annotated as
      Goblet cells, Ciliated Cells, and Ciliated Cells (Figure 3-13). The remaining cells exhibited high expression of
      <em>SFTPB</em> (Figure 3-13), allowing them
      to be broadly classified as alveolar cells. The expression of markers (<em>KRT7, MUC5AC, MUC5B</em> for
      Goblet cells; <em>FOXJ1, CCDC17, TUBB4B</em> for Ciliated cells and <em>SFTPB</em> for Alveolar cells) in Figure 3-14 validated the
      reliability of our annotations.
    </p>

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

    <!--Figure 3-11-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-11.jpg" />
      <b class="image-title">Figure 3-11</b>
    </div>

    <!--Figure 3-12-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-12.jpg" />
      <b class="image-title">Figure 3-12</b>
    </div>

    <!--Figure 3-13-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-13.jpg" />
      <b class="image-title">Figure 3-13</b>
    </div>

    <!--Figure 3-14-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-14.jpg" />
      <b class="image-title">Figure 3-14</b>
    </div>


    <p class="body">
      <b>Step 7.</b> Then, we selecte the SFTPB+ cells in Local Plot &#9316; for further analysis. We create Local Plot
      &#9317; with the pipeline parameter in Table 3-6. The Local Plot &#9317; was identified as nine sub-clussters.
      We used Human - Azimuth Cell Types 2021 (for c_5), CellMarker2024 (for c_7), CellMarker2024(for c_8) as the gene set
      databases and Enrichr as the annotation
      recommendation method. We applied logfoldchange thresholds to filter marker genes: c_5: 4,
      c_7: 4, c_8: 4. The cell type recommendations are shown as Figure 3-15, 3-16. The expression of markers (
      <em>SFTPC, SFTPB, SFTPA1, SFTPD</em> for Pulmonary Alveolar Type II Cells and <em>AGER, CAV1, PDPN, CLIC5</em>
      for Pulmonary Alveolar Type I Cells) in Figure 3-15, 3-16 validated the
      reliability of our annotations. Besides, the expression patterns of clusters c_1 and c_6 have not been reported in
      previous studies. Cluster c_1 exhibited dual positive <em>SMIM22+</em> and <em>SFTPB+</em> (Figure 3-17), while cluster c_6 displayed
      club cell
      related genes <em>SCGB3A1-</em> and <em>SCGB3A2+</em> (Figure 3-17). Therefore, we empirically named these two epithelial cell
      subtypes based on
      their marker genes. The annotation results are shown in Figure 3-18.
    </p>


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


    <!--Figure 3-15-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-15.jpg" />
      <b class="image-title">Figure 3-15</b>
    </div>

    <!--Figure 3-16-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-16.jpg" />
      <b class="image-title">Figure 3-16</b>
    </div>

    <!--Figure 3-17-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-17.jpg" />
      <b class="image-title">Figure 3-17</b>
    </div>

    <!--Figure 3-18-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-18.jpg" />
      <b class="image-title">Figure 3-18</b>
    </div>


    <p class="body">
      <b>Step 8.</b> We select the unannotated cells in Local Plot &#9317; for further analysis. We create Local Plot
      &#9318; with the pipeline parameter in Table 3-7. The Local Plot &#9318; was identified as six sub-clussters. We
      used Human - PanglaoDB
      Augmented 2021 as the gene set database and Enrichr as the annotation
      recommendation method. We applied logfoldchange threshold to filter marker genes: c_5: 1. The recommendation result
      are shown as Figure 3-19. Therefore, c_5 was annotated as basal cells.

      Accorrding to the high expression of BASCs markers (Figure 3-20), c_0 and c_2 can be annotated as BASCs.

      Similarly, c_1 can be annotated as <em>GNB2L1+</em> AT II cells. c_3 can be annotated as <em>HLA-DRB5+SFTPB+</em> ECs amd c_4 can be
      annotated as <em>RACK+SFTPB-ECs</em> (unreported epithelial cell types).
      
      The annotation results are shown as Figure 3-21.

    </p>

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

    <!--Figure 3-19-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-19.jpg" />
      <b class="image-title">Figure 3-19</b>
    </div>

    <!--Figure 3-20-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-20.jpg" />
      <b class="image-title">Figure 3-20</b>
    </div>

    <!--Figure 3-21-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-21.jpg" />
      <b class="image-title">Figure 3-21</b>
    </div>

    <p class="body">
      <b>Step 9.</b> Returning to the Global Plot &#9312;, we next select c_3 to generate a new Local Plot &#9319; (the
      pipeline
      parameters are listed in Table 3–8).
      In Local Plot &#9319;, we identified three sub-clusters. We used Human - PanglaoDB
      Augmented 2021 as the gene set database and Enrichr as the annotation
      recommendation method to annotate each group. We applied
      logfoldchange thresholds to filter marker genes for each group: c_0: 5,
      c_1: 3, c_2: 5. The cell type recommendations are shown as Figure 3-22. Therefore, these sub-clusers were annotated
      as
      Endothelial cells, Mesenchymal cells and Leukocytes (Figure 3-23).
    </p>

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

    <!--Figure 3-22-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-22.jpg" />
      <b class="image-title">Figure 3-22</b>
    </div>

    <!--Figure 3-23-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-23.jpg" />
      <b class="image-title">Figure 3-23</b>
    </div>


    <p class="body">
      <b>Step 10.</b> After all annotation operations, the hierarchical tree is shown as Figure 3-24.

      we first merge the Local &#9314; to the Local Plot &#9313;, and this merge option includes
      "Merge
      with Labels(Annotations)", "Merge Labels(Annotations) with the same name"
      , and "Merge small Labels(Annotations)".


      Next, we merge the Local Plot &#9318; to the Local Plot &#9317; ,and then merge the Local Plot &#9317; to the Local
      Plot &#9316;. These merge options includes "Merge
      with Labels(Annotations)", "Merge Labels(Annotations) with the same name"
      , and "Merge small Labels(Annotations)".

      Thrid, we merge the Local Plot &#9313;, the Local Plot &#9315; and then merge the Local Plot &#9316; to the Global Plot &#9312;. This merge option includes "Merge
      with Labels(Annotations)", "Merge Labels(Annotations) with the same name"
      , and "Merge small Labels(Annotations)".

      At last, we merge the Local Plot &#9319; to the Global Plot &#9312;.
      This merge option includes "Merge with Projection" , "Merge
      with Labels(Annotations)", "Merge Labels(Annotations) with the same name"
      , and "Merge small Labels(Annotations)".

      The merged result are shown as
      Figure 3-25.
    </p>


    <!--Figure 3-24-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-25.jpg" />
      <b class="image-title">Figure 3-24</b>
    </div>

    <!--Figure 3-25-->
    <div class="image-container">
      <img class="image" src="tutorials/img/repo-Figure3-24.jpg" />
      <b class="image-title">Figure 3-25</b>
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