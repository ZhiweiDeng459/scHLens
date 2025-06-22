<template>
        <!--上传Gene Sets对话框-->
        <el-dialog
        title="Add a custom Gene Set"
        :visible.sync="showDialog"
        :show-close="true"
        :lock-scroll="true"
        :modal-append-to-body="false"
        width="400px">
        <div style="
            display:flex;
            width:100%;
            height:100%;
            flex-direction:column;
            align-items:center;
            justify-content: space-between;
        ">
            <el-upload
                drag
                multiple
                ref="upload"
                :limit="1"
                :file-list="fileList"
                :accept="'application/json'"
                :on-change="uploadFileChangeHandler"
                :on-remove="uploadFileRemoveHandler"
                :on-success="UploadFileSuccessHandler"
                :on-error="UploadFileErrorHandler"
                :on-exceed="UploadFileExceedHandler"
                :data="{
                    'name':geneSetName,
                    'type':geneSetType,
                    'JobId':JobId
                }"
                :auto-upload="false"
                :action="uploadPath">
                <i class="el-icon-upload"></i>
                <div class="el-upload__text">Drag file to this,or &nbsp;<em>Click to add</em></div>
            </el-upload>

            <div>
                <div class="el-upload__tip" slot="tip" style="font-family:YaHei;font-size:16px;">Set gene set type and name:</div>
                <div
                    style="
                    display:flex;
                    width: 360px;
                    margin-top:10px;
                    justify-content:space-between
                    ">

                    <el-input
                        placeholder="Name"
                        v-model="geneSetName"
                        style="width:220px">
                        <el-select
                                v-model="geneSetType"
                                placeholder="Select..."
                                size="small"
                                style="width:100px;"
                                slot="prepend"
                                >
                                <el-option v-for="item in geneSetTypeList" :key="item" :label="item" :value="item"> </el-option>
                        </el-select>
                    </el-input>
                    <el-button size="small" type="success" @click="submitUpload">upload to server</el-button>
                </div>
            </div>
        </div>
    </el-dialog>

</template>

<script>
export default {
    name:"UploadGeneSetsPanel",
    props:['gene_sets_info','updateGeneSets'],
    computed:{
        JobId(){
            return this.$store.state.JobId
        },
    },
    data(){
        return {
            /**
             * 关键数据
             */
             showDialog:false,//对话框是否显示
             fileList:[],//文件的List
             geneSetName:'',//基因集名词
             geneSetType:'Human',//基因集种类
             geneSetTypeList:['Human','Mouse'],
             /**
              * gene set上传相关
              */
             uploadPath:'/api/uploadGeneSets',

        }
    },
    methods:{
        /**
         * 外部接口
         */

         open(){//打开对话框
            const self = this;
            self.showDialog = true;

        },

        close(){//关闭对话框
            const self = this;
            self.showDialog = false;
        },


        /**
         * 内部函数
         */

        /**
         * 文件操作handler
         */
        uploadFileChangeHandler(file,fileList){
            this.fileList = fileList
        },
        uploadFileRemoveHandler(file,fileList){
            this.fileList = fileList
        },
        UploadFileSuccessHandler(res,file,fileList){
            this.$message({
                'message':'UPLOAD COMPLETE',
                'type':'success',
                'showClose':true,
            })
            this.$refs.upload.clearFiles();
            this.fileList = []
            this.geneSetType = 'Human'
            this.geneSetName = ''
            this.close();

            //更新gene sets
            this.updateGeneSets();

        },
        UploadFileErrorHandler(err,file,fileList){
            this.$message({
                'message':'UPLOAD FAILED',
                'type':'error',
                'showClose':true,
            })
            this.$refs.upload.clearFiles();
            this.fileList = []
            this.geneSetType = 'Human'
            this.geneSetName = ''

        },
        UploadFileExceedHandler(file,fileList){
            this.$message({
                'message':`Only one gene set file can be uploaded at a time.`,
                'type':'info',
                'showClose':true,
            })
        },
        submitUpload(){
            //check1：是否有文件检查
            if(this.fileList.length == 0){
                this.$message({
                    'message':`Please select the gene set file`,
                    'type':'error',
                    'showClose':true,
                })
                return;
            }


            //check2：是否输入了gene set名称
            if(this.geneSetName == '' || this.geneSetName === null || this.geneSetName === undefined){
                this.$message({
                    'message':"Please input a valid gene set's name",
                    'type':'error',
                    'showClose':true,
                })
                return;
            }
            //check3：重名检查
            for(let type_item of this.gene_sets_info){
                if(type_item['value'] == this.geneSetType){
                    for(let temp_set of type_item['children']){
                        if(temp_set['value'] == this.geneSetName){
                            this.$message({
                                'message':"This gene set name has been used",
                                'type':'error',
                                'showClose':true,
                            })
                            return;
                        }
                    }
                }
            }
            //check4：名称是否包含非法字符检查
            let forbiddenCharas = ['/','+','-','.']
            for(let c of forbiddenCharas){
                if(this.geneSetName.includes(c)){
                    this.$message({
                        'message':`Invalid name: contains forbidden character "${c}"`,
                        'type':'error',
                        'showClose':true,
                    })
                    return;
                }
            }


            //上传
            this.$refs.upload.submit();

        },

    },
    mounted(){
        this.uploadPath = `${window.location.pathname}api/uploadGeneSets`
    }
}
</script>

<style scoped>

</style>