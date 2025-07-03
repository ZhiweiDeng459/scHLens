<template>
  <div class="container">
    <b class="title">Analysis with Views</b>
    <p class="body">
      Upon completion of the pipeline, scHLens generates a variety of
      interactive views that provide a comprehensive profile of the data. Aside
      from the views, there is also a controller bar to control these views.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/system_layout.png" />
      <b class="image-title">Figure 1. The Components' Location in scHLens</b>
    </div>
    <p class="body">
      Figure 1 illustrates the positions of these views and the controller
      within the system: Areas A-E correspond to <b>Cell Projection View</b>,
      <b>Gene Projection View</b>, <b>Distribution of Gene Expression View</b>,
      <b>Marker Gene View</b> and <b>Hierarchy Navigation</b>. Area F is
      <b>controller bar</b>. The following sections will provide a detailed
      introduction to these components and their associated interactions.
    </p>

    <h1 id="t_Views">Views</h1>
    <p class="body">
      In this section, we will introduce the visual representations of views.
    </p>

    <h2 id="t_CellPrjectionView">Cell Projection View</h2>
    <div class="image-container">
      <img class="image" src="tutorials/img/CellProjection.jpg" />
      <b class="image-title">Figure 2. The Cell Projection View</b>
    </div>
    <p class="body">
      Just shown as Figure 2, The Cell Projection View represents the
      distribution of cells in 2-dimensional space using a scatter plot (The
      positions of the cells are obtained from the visualization module within
      the pipeline). Each dot in the plot represents a cell, and the color of
      the dot indicates its group from the clustering module in the pipeline.
      The legend on the right also shows the correspondence between group colors
      encoding and group annotations.
    </p>
    <p class="body">
      Besides, <b>the silhouette score</b>, a computed cluster assessment score,
      will be displayed in the top-left corner of the view.
    </p>
    <p class="body">
      Selected cells will be highlighted with a black border in this view. For
      more information about cell selection, refer to the section "Tutorial ->
      Analysis with Views -> Interaction -> Select Cells"
    </p>

    <h2 id="t_GeneProjectionView">Gene Projection View</h2>
    <div class="image-container">
      <img class="image" src="tutorials/img/GeneExpression.jpg" />
      <b class="image-title">Figure 3. The Gene Projection View</b>
    </div>
    <p class="body">
      The Gene Projection View represents the distribution of specific gene(s)
      expression on a 2-dimensional scatter plot, with points positioned
      identically to those in the Cell Projection View (as shown in Figure 3).
      It consists of two modes: "Scatter" and "Density". Analysts can switch
      mode using the switcher located in the top-right corner of the view (Area
      A in Figure 3).
    </p>
    <ol>
      <li class="body">
        <b>In "Scatter" mode</b>, a gradient color bar represents expression
        values, with each cell dot's color indicating its corresponding value
        along the right color legend.
      </li>
      <li class="body">
        <b>In "Density" mode</b>, the view uses a heatmap to illustrate the
        density of expression values. The right color legend also represents the
        correspondence between the colors and the expressed values. The slider
        (Area B in Figure 3(B)) controls the resolution of the heatmap. A higher
        resolution will make the heatmap display finer expressions, but it may
        also lead to the plot being fragmented.
      </li>
    </ol>
    <p class="body">
      The "Current Gene" component determines the genes to be displayed in this
      view, which will be then explained in detail in the section "Tutorial ->
      Analysis with Views -> Interactions -> Switch Current Genes".
    </p>

    <h2 id="t_DiffView">Distribution of Gene Expression View</h2>
    <div class="image-container">
      <img class="image" src="tutorials/img/diff_gene.jpg" />
      <b class="image-title"
        >Figure 4. The Distribution of Gene Expression View</b
      >
    </div>
    <p class="body">
      The Distribution of Gene Expression View uses a violin plot to indicate
      the expression of specific gene(s) in each group (as shown in Figure 8),
      aiming to compare gene expression differences across groups. These groups
      are obtained from the clustering module in the pipeline. The violin plot
      uses kernel density plot to show the distribution of the gene's expression
      across groups.
    </p>
    <p class="body">
      The plot also offers two modes to assist analysts: "Box" mode and "Strip"
      mode. Analysts can switch mode using the switcher located in the top-right
      corner of the view (Yellow box in Figure 4).
    </p>
    <ol>
      <li class="body">
        <b>In "Box" mode</b>, a white dot represents the median, and the central
        black box indicates the interquartile range. The thin black line
        delineates the remaining distribution.
      </li>
      <li class="body">
        <b>In "Strip" mode</b>, the view visualizes all cells as black dots on
        the kernel density plots.
      </li>
    </ol>
    <p class="body">
      The "Current Gene" component determines the genes to be displayed in this
      view, which will be then explained in detail in the section "Tutorial ->
      Analysis with Views -> Interactions -> Switch Current Genes".
    </p>
    <p class="body">
      If you have selected some cells ("Tutorial -> Analysis with Views ->
      Interaction -> Select Cells"), the "Strip" mode will display only the
      black dots representing the selected cells.
    </p>

    <h2 id="t_MarkerGeneView">Marker Gene View</h2>
    <div class="image-container">
      <img class="image" src="tutorials/img/marker_gene.jpg" />
      <b class="image-title">Figure 5. The Marker Gene View</b>
    </div>
    <p class="body">
      The Marker Gene View uses a dot plot to indicate the marker genes'
      expression in each group (as shown in Figure 5). These marker genes are
      obtained from the DEG Identification module in the pipeline.
    </p>
    <p class="body">
      In the dot plot, the x-axis represents the marker gene list, and the
      y-axis represents the group list. Each dot indicates the expression levels
      of the marker gene in the group. Dot size represents the proportion of
      expressed cells, and the color gradation of dots represents the mean
      expression value.
    </p>
    <p class="body">
      You can also left-click on a marker gene on the x-axis to add it to the
      current genes list.
    </p>

    <h2 id="t_HierarchyNavigation1">Hierarchy Navigation</h2>
    <p class="body">
      We will introduce this part in detail in the chapter "Hierarchy
      Exploration".
    </p>

    <h2 id="t_Controller">Controller Bar</h2>
    <p class="body">
      Area F in Figure 1 is the <b>Controller Bar</b>, which facilitates general
      operations. The Controller Bar displays the current cell count and gene
      count (as shown in Figure 6) and includes functionalities such as Switch
      Current Genes, Edit Group Information, and Refreshing Views. These
      features will be explained in detail in the section "Tutorial -> Analysis
      with Views -> Interactions".
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/count-info.png" />
      <b class="image-title">Figure 6. current cells and genes information</b>
    </div>

    <h1 id="t_Interactions">Interactions</h1>

    <h2 id="t_SelectCells">Select Cells</h2>
    <p class="body">
      Selecting cells is often a prerequisite for other interactions, such as
      Export Cells. The selected cells are highlighted with black borders in the
      Cell Projection View. scHLens provides multiple methods for selecting
      cells:
    </p>
    <ol>
      <li class="body">
        <b>Lasso Selection</b>: In the Cell Projection View, you can drag the
        mouse to create a lasso to select cells of interest.
      </li>
      <li class="body">
        <b>Select by Group</b>: In the Cell Projection View, click on a group in
        the legend to select all cells in that group.
      </li>
      <li class="body">
        <b>Select by Gene Filter</b>: scHLens provides a Gene Filter to assist
        users to select cells. This filter can select cells that meet specific
        gene expression conditions. Shown as Figure 7, Users can click blue
        "Filter Cells" button in the left side of the system to open Gene
        Filter. Users can type a gene name or part of it in the input box and
        select a matching gene to add to the filter. And then they can drag the
        slider to set the range of genes as conditions. Finally, they can click
        the blue button "Filter" to select cells that meet specific gene
        expression conditions.
      </li>
      <li class="body">
        <b>Select All</b>: In the context menu of the Cell Projection View, you
        can choose "Select All" to select all cells.
      </li>
      <li class="body">
        <b>Inverse Selection</b>: In the context menu of the Cell Projection
        View, you can choose "Inverse Selection" to invert the selection of
        cells (selecting unselected cells and deselecting previously selected
        ones).
      </li>
      <li class="body">
        <b>Clear Selection</b>: In the context menu of the Cell Projection View,
        you can choose "Clear Selection" to deselect all cells.
      </li>
    </ol>
    <div class="image-container">
      <img class="image" src="tutorials/img/Filter_cells.png" />
      <b class="image-title">Figure 7. Cells Filter Tool</b>
    </div>

    <h2 id="t_DeleteCells">Delete Cells</h2>
    <p class="body">
      To prevent some unimportant cells from disturbing the analysis, Deleting
      some cells is essential sometimes. In the Cell Projection View, you can
      click "Delete Chosen Elements" in the context menu after selecting them.
    </p>

    <h2 id="t_ExportCells">Export Cells</h2>
    <p class="body">
      When selecting some important cells, saving them as a single dataset for
      analysis in other scenarios is a good choice. To do this, select the cells
      and click the "Save as Dataset" option in the context menu in the Cell
      Projection View. Then, you should name the new dataset in the popup
      dialog, and the dataset will be exported in h5ad format.
    </p>

    <h2 id="t_MergeDuplicateLabels">Merge Duplicate Labels(Annotations)</h2>
    <p class="body">
      During annotation (refer to the section "Tutorial -> Analysis with Views
      -> Interaction -> Annotation "), different groups may be labeled as the
      same cell type. In this case, merging groups with identical labels into
      one group is essential. In the Cell Projection View, you can click the
      "Merge Duplicate Labels" option in the context menu, and groups with the
      same labels will be merged into the same group.
    </p>

    <h2 id="t_Recovery">Recovery</h2>
    <p class="body">
      During hierarchical exploration(See Section "Hierarchical Exploration" for
      details), the position and labels(annotations) of cells may change. If you
      want to restore the data to its initial state, you just need to click the
      “Restore the Labels” option or the “Restore the Projection” option in the
      context menu of the Cell Projection View."
    </p>

    <h2 id="t_SwitchGenes">Switch Current Genes</h2>
    <p class="body">
      Gene Projection View and Distribution of Gene Expression View display
      information about one or several specific genes (<b
        >We call them "Current Gene"</b
      >). These genes are shown in the Controller Bar (Figure 8) and can be
      changed.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/CurrentGenes.png" />
      <b class="image-title">Figure 8. Current Genes</b>
    </div>
    <p class="body">
      Shown in Figure 8, the "Current Gene" component displays the list of
      currently selected genes. Clicking the cross icon to the right of a gene
      name removes the gene from the list. The red "Reset" button can clear the
      list and resets it to the default gene.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/SelectGenes.png" />
      <b class="image-title">Figure 9. Select Genes</b>
    </div>
    <p class="body">
      scHLens supports multiple ways to select genes for addition to the current
      gene list. As shown in Figure 9(A), users can type a gene name or part of
      it in the input box in Controller Bar and select a matching gene to add to
      the current gene list. If users want to add multiple genes at once, they
      can click the button highlighted in the red box in Figure 9(A), which
      opens the "Multiple Genes Input" component in Figure 9(B). In this
      component, users can efficiently add multiple genes by entering them
      separated by semicolons, spaces, or line breaks.
    </p>

    <h2 id="t_EditGroup">Edit group color & Export group markers</h2>
    <p class="body">
      As shown in Figure 10, the blue button "Group" in Controller Bar will open
      a dialog with group information(<b>Group Info dialog</b>). These groups
      are generated by the clustering algorithm selected by analysts in the
      Cluster module of the pipeline. The Group Info dialog displays information
      such as the group's ID, color, size, and annotation(label). The Group Info
      dialog also supports the following interactions:
    </p>
    <ol>
      <li class="body">
        <b>Edit group color</b>: Users can click the group color mark(Such as
        Area A in Figure 10) to personalize its color.
      </li>
      <li class="body">
        <b>Export group markers</b>: Users can click the group's blue download
        button (such as Area B in Figure 10) to download its marker genes
        generated by the DEG Identification module in the pipeline. The
        downloaded files includes two types: global markers and local markers.
        Global markers refer to the marker genes calculated by placing the
        current group at the global root node. Local markers refer to the marker
        genes calculated within the pipeline of the current node.(Regarding the
        global node and current node, see Section "Hierarchical Exploration" for
        details.)
      </li>
    </ol>
    <div class="image-container">
      <img class="image" src="tutorials/img/Group_info.png" />
      <b class="image-title">Figure 10. The Group Info Dialog</b>
    </div>

    <h2 id="tAnnotation">Annotation</h2>
    <p class="body">
      Cell type annotation is one of the core tasks of scHLens. scHLens supports
      annotation based on the <b>Group Info dialog</b> (see "Tutorial ->
      Analysis with Views -> Interactions -> Edit group color & Export group
      markers"). Users can annotate manually or based on automated
      recommendations.
    </p>
    <p class="body">
      For <b>manual annotation</b>, users only need to edit the Annotation text
      box in <b>Group Info dialog</b>. For
      <b>annotation based on automated recommendations</b>, the process is shown
      in Figure 11(A):
    </p>
    <ol>
      <li class="body">
        Click the text box of the group in Group Info dialog and then Pop the Recommendation dialog.
      </li>
      <li class="body">
        Choose the recommendation method and the gene set library
      </li>
      <li class="body">
        Click the "Filter" button and configure the recommendation parameters:
        The "Gene Filter" helps filter marker genes used for cell
        recommendations. For example, the "logfoldchanges" option filters out
        marker genes with low log fold change values. The "Result Filter" helps
        users filter the cell recommendation results.
      </li>
      <li class="body">
        Click green button "Run" and get the recommendation results.
      </li>
      <li class="body">
        Select an appropriate cell type in the result table and click it to
        automatically replace the group's annotation.
      </li>
    </ol>
    <p class="body">
      The recommended results are displayed both in tabular form(Figure 11(A))
      and as a histogram ( Figure 11(B)).
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/Annotation.jpg" />
      <b class="image-title">Figure 11. Annotation</b>
    </div>
    <p class="body">
      In addition to using the built-in Gene Sets, you can also upload your own
      Gene Set. As shown in Figure 12, click the blue "+" button and upload a
      JSON-formatted Gene Set file in the pop-up window. The Gene Set format
      follows the example of
      <a
        href="http://www.csuligroup.com:32768/resource/sampleGeneSets.json"
        download="Sample.json"
        >Sample.json</a
      >.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/uploadGeneSet.png" />
      <b class="image-title">Figure 12. Upload your Gene Set</b>
    </div>
    <h2 id="t_SaveImages">Save the images</h2>
    <p class="body">
      Each view, except for the Hierarchy Navigation, can be saved as an
      SVG file. To do this, you need to click the "Save this image" option in
      the context menu of these views.
    </p>
    <h2 id="t_Refreshing">Refreshing Views</h2>
    <p class="body">
      The right button in Controller Bar, marked with a recycling icon, is used
      to refresh all views with a single left-click.
    </p>
  </div>
</template>

<script>
export default {
  name: "TAnalysiswithViews",
  methods: {},
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