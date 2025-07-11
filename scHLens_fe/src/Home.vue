<template>
  <div class="home-container">
    <div class="home-row1">
      <div class="content-box content1">
        <div class="content-box-title">
          <b class="content-box-title-text">Visit Visualization</b>
        </div>
        <!-- <a href="https://info.flagcounter.com/Wrtg"
          ><img
            style="margin: 20px"
            src="https://s01.flagcounter.com/map/Wrtg/size_m/txt_000000/border_CCCCCC/pageviews_0/viewers_0/flags_0/"
            alt="Flag Counter"
            border="0"
        /></a> -->
        <!-- <a href='https://clustrmaps.com/site/1c56d'  title='Visit tracker'>
          <img src='//clustrmaps.com/map_v2.png?cl=ffffff&w=550&t=tt&d=2KYfHIZwB-VTu7d5cWrNgBIQJMtFcIzN8wpfedn98qM&co=2d78ad&ct=ffffff'/>
        </a> -->
        <img
          src="//clustrmaps.com/map_v2.png?cl=ffffff&w=550&t=tt&d=2KYfHIZwB-VTu7d5cWrNgBIQJMtFcIzN8wpfedn98qM&co=2d78ad&ct=ffffff"
        />
      </div>
      <div class="content-box content2">
        <div class="content-box-title">
          <b class="content-box-title-text">What is scHLens?</b>
        </div>
        <div class="content-box-content">
          <a class="content-box-content-text">
            <b>scHLens</b> is a web server for hierarchically and interactively
            exploring cellular heterogeneity. scHLens focused on the cell type
            annotation of scRNA-seq. The main novelty of scHLens includes two
            aspects. Firstly, for each layer, scHLens provided a customized
            pipeline for gene selection, cell visualization, clustering,
            differential expressed gene analysis and cell type annotation. Then,
            scHLens provided a hierarchical analytical mode to explore cellular
            heterogeneity, and proposed a novel hierarchically visualization
            method to integrate different levels of cells subgroups. This web
            platform allows users to easily select, export, delete and reanalyze
            specific cell groups of their interest. We believe that this
            hierarchical explo-ration of single-cell RNA-seq data will become a
            new paradigm for the analysis of cellular heterogeneity
          </a>
        </div>
      </div>
    </div>
    <div class="home-row2">
      <div class="content-box content3">
        <div class="content-box-title">
          <b class="content-box-title-text">Getting Started with scHLens</b>
        </div>
        <div class="content-box-content" style="flex-direction: column">
          <a class="content-box-content-text">
            <a><a style="font-size: 15px">▶&nbsp;</a> In scHLens, You can:</a>
            <ul>
              <li>
                To learn how to use scHLens, click the "Tutorial" button in the
                top menu to access the tutorials.
              </li>
              <li>
                To start using scHLens directly, click the "System" button in
                the top menu to enter the system.
              </li>
            </ul>
          </a>
          <a class="content-box-content-text">
            <a style="font-size: 15px">▶&nbsp;</a> Beside, scHLens have also
            intergrated the three cases in the manuscript. For more details
            about how to use the three cases, please refer to: "Tutorials ->
            Enter a Job -> Create a New Job & Load an Existing Job & Use Samples
            -> Use Samples"
          </a>
          <a class="content-box-content-text"
            ><a style="font-size: 15px">▶&nbsp;&nbsp;</a>There is a
            demonstration video available for scHLen (However, this video only briefly 
            demonstrates the main workflow of the system rather than all
            functions, in order to avoid making it too lengthy. For more
            details, please refer to the tutorials.) :</a
          >
          <video
            controls
            style="
              margin-top: 20px;
              width: 60%;
              height: auto;
              display: block;
              margin-left: auto;
              margin-right: auto;
            "
            :src="'http://www.csuligroup.com:32768/resource/Demovideo.mp4'"
          >
            Your browser does not support this video. Please use the latest
            version of Chrome or Edge.
          </video>
        </div>
      </div>
    </div>

    <div class="home-row2">
      <div class="content-box content3">
        <div class="content-box-title">
          <b class="content-box-title-text">The framework of scHLens</b>
        </div>
        <div class="content-box-content">
          <img
            src="pics/framework.jpg"
            style="width: 100%; height: 100%; margin: 20px"
          />
          <!-- <a class="content-box-content-text">
            Genome-wide association studies (GWAS) have identified thousands of
            noncoding loci that are associated with human diseases and complex
            traits, each of which could reveal insights into the mechanisms of
            disease (5). Many of the underlying causal variants may affect
            enhancers, but we lack accurate maps of enhancers and their target
            genes to interpret such variants (6-7). Thus, the emerging
            computational methods facilitate the identification of E-P
            interactions, as well as the interpretation of pathological
            variants.
          </a> -->
        </div>
      </div>
    </div>
    <!-- <div class="home-row3">
      <div class="content-box content3">
        <div class="content-box-title">
          <b class="content-box-title-text"
            >Why E-P interactions are important?</b
          >
        </div>
      </div>
    </div> -->
    <div class="home-end"></div>
  </div>
</template>

<script>
export default {
  name: "Home",
  components: {},
  computed: {
    curData() {
      return this.$store.state.curData;
    },
    dataList() {
      return this.$store.state.dataList;
    },
    pipelineEntityArray() {
      return this.$store.state.pipelineEntityArray;
    },
    JobId() {
      return this.$store.state.JobId;
    },
  },
  mounted() {},
};
</script>

<style scoped lang="less">
.home-container {
  background-color: #f5f5f5;
  display: flex;
  flex: 1 1 0;
  flex-direction: column;
  align-items: center;
  overflow-x: hidden;
  overflow-y: auto;
  width: 100%;
  .content-box {
    background-color: white;
    margin: 10px;
    border: 2px solid rgb(200, 200, 200);
    border-radius: 4px;
    box-shadow: 0px 0px 9px lightgray;
    display: flex;
    flex-direction: column;
    align-items: center;
    .content-box-title {
      margin: 20px 10px 10px 10px;
      .content-box-title-text {
        font-size: 30px;
        border-bottom: 4px solid black;
        padding-bottom: 6px;
      }
    }
    .content-box-content {
      flex: 1 1 0;
      display: flex;
      align-items: flex-start;
      justify-content: center;
      padding: 20px 60px;
      .content-box-content-text {
        font-size: 18px;
        line-height: 30px;
        margin-bottom: 10px;
      }
    }
  }
  .home-row1 {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    margin-top: 30px;
    width: 80%;
    .content1 {
    }

    .content2 {
      width: 62%;
    }
  }
  .home-row2 {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    margin-top: 30px;
    width: 80%;
    .content3 {
      width: 100%;
    }
  }
  .home-row3 {
    margin-top: 30px;
    width: 80%;
  }
  .home-end {
  }
}
</style>