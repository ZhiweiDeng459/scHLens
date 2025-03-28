<template>
  <div class="container">
    <b class="title">Select a Dataset</b>
    <p class="body">
      After entering the job, you need to select a dataset as the target of the
      analysis before setting up the pipeline. This can be done by choosing from
      the provided sample datasets in the dropdown menu (Area A in Figure 1) at
      the top-left of scHLens. You can also add your own dataset by clicking the
      'plus' button (Area B in Figure 1). Area C in Figure 1 shows the size of
      the selected dataset.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/Dataset.png" />
      <b class="image-title">Figure 1. Select The Dataset</b>
    </div>

    <h1 id="t_SelectSample">Select Sample Dataset</h1>
    <div class="image-container">
      <img class="image" src="tutorials/img/Select_sample_dataset.png" />
      <b class="image-title">Figure 2. Select Sample Dataset</b>
    </div>
    <p class="body">
      scHLens provides three sample datasets for analysis, including "3k-PBMCs",
      "LUAD", "muraro" (as shown in Figure 2).
    </p>

    <h1 id="t_AddDataset">Add your own Dataset</h1>
    <div class="image-container">
      <img class="image" src="tutorials/img/Add_custom_dataset.jpg" />
      <b class="image-title">Figure 3. Add Custom Dataset</b>
    </div>
    <p class="body">
      scHLens supports datasets in various formats, including csv, MEX, h5ad,
      and h5df. To add your own dataset, you should follow these steps:
    </p>
    <ol>
      <li class="body">
        Select the dataset format for upload (Area C in Figure 3).
      </li>
      <li class="body">
        Read the format-specific upload requirements (Area A in Figure 3).
      </li>
      <li class="body">
        After preparing the dataset for upload, analysts can either drag the
        dataset file to the area B in Figure 3, or click "Click to add" in the
        area B in Figure 3 and navigate to the location of the local dataset
        file to complete the upload process.
      </li>
      <li class="body">
        Name this dataset in the textbox (Area C in Figure 3).
      </li>
      <li class="body">
        Click "upload to server" to complete the dataset upload (Area D in
        Figure 3).
      </li>
      <li class="body">
        You can find this dataset in the dataset dropdown menu.
      </li>
    </ol>
    <a class="body">
      <b class="body">Note:</b>&nbsp;&nbsp;The "Stratified" method in the
      Downsampling module of the pipeline requires a dataset with labels. If you
      want to upload a dataset with labels, you should upload it in h5ad format
      and set the labels in AnnData.obs['label'].
    </a>
  </div>
</template>

<script>
export default {
  name: "TselectDataset",
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