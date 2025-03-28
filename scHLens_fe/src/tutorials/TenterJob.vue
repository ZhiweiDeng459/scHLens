<template>
  <div class="container">
    <b class="title">Enter a Job</b>
    <p class="body">
      Upon reaching the "System" page, a dialog will pop up (Area A in Figure 1)
      to notify you should start analysis within a job. You can create a new job
      or load an existing job.
    </p>
    <div class="image-container">
      <img class="image" src="tutorials/img/Job_enter.png" />
      <b class="image-title">Figure 1. The Job Dialog</b>
    </div>

    <h1 id="t_WhatisJob">What is Job?</h1>
    <p class="body">
      A job is a personal workspace that records the user's analysis process and
      results for a scRNA-seq dataset. The job mechanism allows users to
      reproduce their previous work the next time they open the system, even
      after closing it.
    </p>

    <h1 id="t_WhatisJobID">What is Job ID?</h1>
    <p class="body">
      Each job is assigned a unique "Job ID", which acts as its identifier. This
      ID is crucial for loading previous jobs.
    </p>

    <h1 id="t_CreateLoadJob">
      Create a New Job & Load an Existing Job & Use Samples
    </h1>
    <p class="body">
      <b>Create a New Job</b>: To begin a new analysis, click the green "Create
      new Job" button (Area C in Figure 1). This action will initiate a new job,
      ready to record the user's analysis process and results.
    </p>
    <p class="body">
      <b>Load an Existing Job</b>: To reproduce a previous work, enter the
      corresponding Job ID into the designated text box (Area B in Figure 1) and
      click the blue "Load Exist Job" button (Area D in Figure 1) . Then this
      job will be loaded, allowing for continuation from the last auto-saved
      state.
    </p>
    <p class="body">
      <b>Use Samples</b>: scHLens provides three Sample Jobs to illustrate the
      cases in manuscript. You can click the 'Use Samples' button (Bottom-Right
      corner in Figure 1) to select a sample to load.
    </p>
    <h1 id="t_ExportImportJob">
      Export the Job & Import the Job
    </h1>
    <div class="image-container">
      <img class="image" src="tutorials/img/Job.png" />
      <b class="image-title">Figure 2. Export / Delete job</b>
    </div>
    <p class="body">
      If you want to migrate a job to another deployment version of scHLens or
      save a job's current state for future rollback, you can export the current
      state of the job as an importable local file. To export the job, click the
      blue 'Export Job' button (Area A in Figure 2), and you will receive a
      compressed job file. To import the local file into a different system,
      click the blue 'Import Saved Job' button (Area E in Figure 1), select the
      file, and you will be able to access this job in the new system via the
      same Job ID. If a job with the same name already exists in the new system,
      the imported job will overwrite the original one.
    </p>

    <h1 id="t_DeleteJob">Delete Job</h1>
    <p class="body">
      To ensure data security, you can delete your job after completing the
      analysis. You just need to click the red 'Delete Job' button (Area B in
      Figure 2).
    </p>
  </div>
</template>

<script>
export default {
  name: "TenterJob",
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