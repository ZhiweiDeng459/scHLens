MySCTransform <- function(matrix,gene_name,cell_name,topGenes){

    library(Seurat)
    library(glmGamPoi)
    library(sctransform)

    rownames(matrix) = cell_name
    colnames(matrix) = gene_name

    rna <- CreateSeuratObject(counts = t(matrix))

    rna <- SCTransform(rna,variable.features.n=topGenes,verbose = FALSE,seed.use=1448145)
    
    result <- as.data.frame(rna@assays$SCT@scale.data)

    return(result)

}

library(Matrix)
library(jsonlite)

# 保存的路径

args <- commandArgs(trailingOnly = TRUE)
JobId <- args[1]
save_path = file.path(".", "job_attachment", JobId, "scTransform")
json_path <- file.path(save_path, "params.json")

# 读取矩阵
mat <- readMM(file.path(save_path, "matrix.mtx"))
mat <- as(mat, "CsparseMatrix")
  
# 读取行名、列名
rownames <- read.table(file.path(save_path, "rows.tsv"), stringsAsFactors = FALSE)[,1]
colnames <- read.table(file.path(save_path, "cols.tsv"), stringsAsFactors = FALSE)[,1]

# 读取其他参数
other_params <-fromJSON(json_path)

# 如果你的 multik 函数定义在 R 里，可以直接调用
result <- MySCTransform(matrix = mat, gene_name = colnames, cell_name = rownames, topGenes = other_params$topGenes)
write.csv(result, file = file.path(save_path, "result.csv"), row.names = TRUE)
