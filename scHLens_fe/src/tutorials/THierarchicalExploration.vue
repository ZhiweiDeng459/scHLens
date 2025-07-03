<template>
  <div class="container">
    <b class="title">Hierarchical Exploration</b>
    <p class="body">
      Because of the hierarchical and complex nature of cell types, a
      comprehensive analysis of scRNA-seq data is challenging to achieve in a
      single analytical run. Therefore, employing a hierarchical and interactive
      exploration way is more appropriate.
    </p>
    <p class="body">
      The workflow of Hierarchical Exploration is simple:
      <b>Annotate groups in the Global Plot</b>, <b>Create Local Plot</b>,
      <b>Find good pattern and Annotate groups in Local Plot</b> and then
      <b>Merge the Local Results</b>. For the process of "Annotate groups in
      Global Plot" and "Annotate groups in Local Plot", please refer to "Tutorial -> Analysis with Views ->
      Interactions -> Annotation." We will describe the remaining
      steps in detail below.
    </p>
    <p class="body">
      <b>Global Plot</b> refers to the initial views created after the first
      pipeline, which includes the largest number of cells. <b>Local Plot</b>
      refers to new views generated after selecting interesting population of
      cells from the Global Plot and running a new pipeline. The Local Plot
      enables the analysis to focus on specific subsets of cells, revealing new
      patterns, such as cell subpopulations, within the new subspace.
      Additionally, cells can be selected from a Local Plot to create subsequent
      Local Plots. This creation process forms a tree-like structure, which is
      represented using Hierarchy Navigation. Detailed explanations can be found
      in the section "Tutorial -> Hierarchical Exploration -> Hierarchical Navigation."
    </p>

    <h1 id="t_CrateLocalPlot">Create Local Plot & Find Good Pattern</h1>
    <p class="body">
      Analysts can select specific cells (as mentioned in "Tutorial -> Analysis with Views ->
      Interactions -> Select Cells") to start a new pipeline and finally create a
      Local Plot.
    </p>
    <p class="body">
      To do this, select the cells, click on the "Create Local Projection"
      option in the context menu of the Cell Projection View (as shown in Figure
      1), and configure the local pipeline as described in the "Start an
      Analysis Pipeline" chapter. By repeatedly performing the "Create Local
      Plot" operation with different pipeline settings, continue until you find
      a good pattern.
    </p>
    <p class="body">
      You can create a new Local Plot not only from the Global Plot but also
      from an existing Local Plot.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/Local_create.jpg" />
      <b class="image-title">Figure 1. Create Local Projection</b>
    </div>

    <h1 id="t_HierarchyNavigation2">Hierarchy Navigation</h1>
    <p class="body">
      The Hierarchy Navigation visualizes the relationship between the Global
      Plot and Local Plots using a tree structure. As shown in Figure 2, the
      root node represents the Global Plot, while the other nodes represent
      Local Plots. Child nodes are created based on the "Create Local Plot"
      operation performed on their parent nodes (see "Tutorial -> Hierarchical Exploration ->
      Create Local Plot & Find Good Pattern" for details). The node
      corresponding to the currently displayed plot in the system interface is
      highlighted in the Hierarchy Navigation, while other undisplayed nodes are
      shown in gray.
    </p>
    <p class="body">
      Hierarchy Navigation supports the following interactions:
    </p>
    <ol>
      <li class="body">
        <b>Overview</b>: Hover over a node to display its Cell Projection View
        as an overview.
      </li>
      <li class="body">
        <b>Switch</b>: To switch the node displayed in the System Interface,
        simply right-click on the corresponding node.
      </li>
      <li class="body">
        <b>Delete</b>: Click the "Cross" button in the top-right corner of the
        node to delete it.
      </li>
    </ol>
    <div class="image-container">
      <img class="image" src="tutorials/img/Hierarchy_Navigation.jpg" />
      <b class="image-title">Figure 2. Hierarchy Navigation</b>
    </div>

    <h1 id="t_DiscoverCorrespondence">
      Discover the correspondence of cells across different nodes
    </h1>
    <p class="body">
      If you want to know the correspondence of cells across different nodes,
      you should choose the "Discover Correspondence" option in the context menu
      of the Cell Projection View (as shown in Figure 3). The Correspondence
      Explorer will pop up after you click this option.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/Correspondence_discover.jpg" />
      <b class="image-title">Figure 3. Discover Correspondence</b>
    </div>
    <p class="body">
      In the Correspondence Explorer(Shown as Figure 4), the left side displays
      a tree identical to the Hierarchy Navigation, while the right side
      provides an overview of the correspondence. You can left-click to select
      an ancestor node of the current node on the tree. The overview on the
      right will then highlight the corresponding positions of the current
      node's cells in the selected node using yellow. If you have previously
      selected a subset of cells, the yellow highlights will only mark the
      selected cells, not all cells in the current node.
    </p>
    <p class="body">
      Notably, if the selected node is not an ancestor of the current node, the
      correspondence cannot be displayed, as their cells do not have an
      inclusion relationship.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/Correspondence_interface.jpg" />
      <b class="image-title">Figure 4. Explorer Interface</b>
    </div>

    <h1 id="t_Merge">Merge the Local Results</h1>
    <p class="body">
      At the end of Hierarchical Exploration, we typically want to merge the
      annotations and good visual patterns from the Local Plot into their
      ancestor nodes (usually the Global Plot). In the Hierarchy Navigation, we
      can left-click to select local nodes as merge-data providers for the merge
      and their ancestor nodes as merge-data receivers. The provider is
      highlighted with a green border (Figure 5(B)), offering annotations and
      good visual patterns; the receiver is highlighted with a blue border
      (Figure 5(A)), merging the provided data into itself.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/merge.jpg" />
      <b class="image-title">Figure 5. Merge Local Result</b>
    </div>
    <p class="body">
      After selecting the nodes, we use the "Merge" option in the context menu
      of the Hierarchical Navigation and set the merge options in the pop-up
      dialog (shown in Figure 6). The "Merge with Projection" option determines
      whether to include the good visual patterns from the Local Plot into the
      merged data. The "Merge with Labels (Annotations)" option determines
      whether to include the annotations from the Local Plot in the merged data.
      The "Merge Labels (Annotations) with the Same Name" option determines
      whether to combine labels with the same name. The "Merge Small Labels
      (Annotations)" option determines whether to merge groups with too few
      cells into their most similar group to avoid issues in calculating marker
      genes.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/merge-options.png" />
      <b class="image-title">Figure 6. Merge Options</b>
    </div>
  </div>
</template>

<script>
export default {
  name: "THierarchicalExploration",
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