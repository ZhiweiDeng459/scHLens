<template>
    <div class="filter-outer-container">
        <div>
            <b class="header-title">Gene Filter</b>
            <HR style="margin:1px 0px"></HR>
        </div>
        <div class="container">
            <div class="control-box">
                <el-autocomplete
                    v-model="state"
                    :fetch-suggestions="querySearch"
                    placeholder="Gene Search"
                    @select="handleSelect"
                    suffix-icon="el-icon-search"
                    size="mini"
                    class="gene-search-box">
                </el-autocomplete>
                <div class="control-option-button" ref="controlOptionButtonInFilter">
                    <!-- <el-popover
                        placement="bottom"
                        trigger="click"
                        v-model="recommendGenesVis"
                        class="recommendGenePopoverInGF">
                            <el-table :data="recommendGenes"  @row-click="selectRecommendGene" style="width:100%" height="300px">
                                <el-table-column class-name="recommend-item" width="200" property="gene" label="Recommended Genes" align="center" style="cursor: pointer;"></el-table-column>
                            </el-table>
                            <el-button  slot="reference" size="mini" type="warning" @click="getRecommendGene" style="width:100%">Recomend</el-button>
                    </el-popover> -->
                    <el-button size="mini" type="danger" @click="clearConditions" class="clearButton">Clear</el-button>
                    <el-button size="mini" type="primary" @click="filterCell" class="filterButton">Filter</el-button>
                </div>
            </div>
            <div class="filter-tip-container">
                <p class="text">{{ filterTip }}</p>
            </div>
            <template v-for="(item, index) in chosenGenes">
                <div :key="item.gene" class="geneItemContainer">
                    <el-checkbox v-model="item.chosen" class="checkbox"></el-checkbox>
                    <span class="text">{{ item.gene }}</span>
                    <el-slider v-model="item.curRange" :min="item.range[0]" :max="item.range[1]" class="slider" range></el-slider>
                    <el-button type="danger" icon="el-icon-delete" size="mini" class="button" @click="deleteSelf(index)" circle></el-button>
                </div>
            </template>
        </div>
    </div>
</template>

<script>
import Vue from "vue";
import { Autocomplete, Checkbox, Button, Slider, Radio ,Popover,Scrollbar} from "element-ui";
import { requestCandidateGeneList, requestGeneValueRange, requestFilteredCellList , recommendGene} from "@/utils/interface.js";

Vue.component(Autocomplete.name, Autocomplete);
Vue.component(Checkbox.name, Checkbox);
Vue.component(Button.name, Button);
Vue.component(Slider.name, Slider);
Vue.component(Radio.name, Radio);
Vue.component(Popover.name,Popover);

export default {
    name: "GeneFilter",
    data() {
        return {
            state: "",
            chosenSet: new Set(),
            chosenGenes: [],
            filterTip: "Search genes to filter cell...",
            recommendGenes:[],
            recommendGenesVis:false,
        };
    },
    computed: {
        cellData() {
            return this.$store.state.curData.cellData;
        },
        dataset(){
            return this.$store.state.curData.paramsObj['dataset'];
        },
        curData(){
            return this.$store.state.curData;
        },
        JobId(){
            return this.$store.state.JobId
        }
    },
    methods: {
        querySearch(queryString, cb) {
            // requestCandidateGeneList(this.JobId,this.curData.ViewId,"^" + queryString)
            //     .then((response) => {
            //         const data = response.data;
            //         cb(
            //             data.map((item) => {
            //                 return { value: item };
            //             })
            //         );
            //     })
            //     .catch((err) => {
            //         console.log(err);
            //     });
            let matched_genes = this.curData.queryGenes.filter(str => str.toLowerCase().startsWith(queryString.toLowerCase()))
            let sorted_matched_genes = matched_genes.sort((a,b)=>a.length - b.length)
            cb(sorted_matched_genes.map(item=>{
                return { value: item };
            }))
        },
        handleSelect(item) {
            if (!this.chosenSet.has(item.value)) {
                this.chosenSet.add(item.value);
                requestGeneValueRange(this.JobId,this.curData.ViewId,item.value)
                    .then((response) => {
                        const range = response.data;
                        range[0] = Math.floor(range[0]);
                        range[1] = Math.ceil(range[1]);
                        this.chosenGenes.push({
                            gene: item.value,
                            range: range,
                            chosen: true,
                            curRange: range,
                        });
                    })
                    .catch((err) => {
                        console.log(err);
                    });
            }
        },
        deleteSelf(index) {
            const geneName = this.chosenGenes.splice(index, 1)[0].gene;
            this.chosenSet.delete(geneName);
        },
        filterCell() {
            requestFilteredCellList(
                this.JobId,
                this.curData.ViewId,
                this.chosenGenes.map((item) => item.gene),
                this.chosenGenes.map((item) => item.curRange)
            )
            .then((response) => {
                const cells = new Set(response.data);
                const chosenCell = this.cellData.filter((item) => cells.has(item.id));
                this.filterTip = `Filtered ${chosenCell.length} cell${chosenCell.length > 1 ? "s" : ""}`;
                this.$store.commit(
                    `updateChosenData`,
                    chosenCell.map((item) => item.id)
                );
            })
            .catch((err) => {
                console.log(err);
            });
        },
        clearConditions(){//清空基因表达条件
            this.chosenGenes.splice(0, this.chosenGenes.length);
            this.chosenSet.clear();

        }
        // getRecommendGene(){
        //     this.recommendGenes.length = 0;
        //     recommendGene(this.JobId,this.curData.ViewId,this.$store.state.recommendMode)
        //     .then((response) => {
        //         for(let item of response.data){
        //             this.recommendGenes.push({
        //                 'gene':item
        //             })
        //         }
        //     })
        //     .catch((err) => {
        //         console.log(err);
        //     });
        // },
        // selectRecommendGene(row){
        //     this.state = row.gene;
        //     this.handleSelect({'value':row.gene})
        //     this.recommendGenesVis = false;
        //     // this.$refs.controlOptionButtonInFilter.click(); //因为el-table没有blur函数，只能用这种粪的方式替代。
        // }
    },
};
</script>

<style scoped lang="less">

.filter-outer-container{
    background-color: white;
    .header-title{
        font-size:20px;
    }
    .container {
        display: flex;
        flex-direction: column;
        align-items: center;
        vertical-align: middle;
        .control-box {
            width: 100%;
            .gene-search-box {
                margin: 10px 10px 0 10px;
                width: calc(100% - 20px);
            }
            .control-option-button{
                margin:5px 10px 2px 10px;
                width: calc(100% - 20px);
                display: flex;
                justify-content:space-between;
                .recommendGenePopoverInGF{
                    width: calc(50% - 5px);
                }
                .filterButton{
                    width: calc(50% - 5px);
                }
                .clearButton{
                    width: calc(50% - 5px);
                }
            }
        }
        .filter-tip-container {
            display: flex;
            width: calc(100% - 12px);
            align-self: start;
            padding: 0 0 5px 12px;
            margin: 7px 0 0 0;
            font-size: 14px;
            color: rgba(150, 150, 150, 0.8);
            border-bottom: 1px solid rgba(190, 190, 190, 0.7);
            .text {
                flex: 1;
                display: inline-block;
                margin: 0;
            }
            .mode-radio-container {
                width: 100px;
                .mode-radio {
                    margin: 0 0 0 10px;
                }
            }
        }
        .geneItemContainer {
            padding-left: 10px;
            width: calc(100% - 10px);
            align-self: start;
            display: flex;
            align-items: center;
            .checkbox {
                margin: 7px 5px 0 0;
            }
            .text {
                flex: 1;
                margin: 8px 0 0 0;
            }
            .slider {
                display: inline-block;
                width: 100px;
                margin: 5px 15px 0 0;
                /deep/ .el-slider__button {
                    width: 12px;
                    height: 12px;
                }
            }
            .button {
                margin: 3px 5px 0 auto;
            }
        }
    }



}

/deep/ .recommend-item{
    cursor: pointer;
}


</style>
