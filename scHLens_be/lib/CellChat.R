CellChat <- function(matrix,cellName,geneName,group,DatabaseType){

    library(Seurat)
    # library(patchwork)
    library(CellChat)
    library(dplyr)

    rownames(matrix) <- cellName
    colnames(matrix) <- geneName
    matrix <- t(matrix)

    meta = data.frame(labels= unlist(group),row.names = unlist(cellName))

    cellchat<-createCellChat(object=matrix,meta = meta,group.by='labels')
    groupSize<-as.numeric(table(cellchat@idents))


    if(DatabaseType=='human'){
        CellChatDB<-CellChatDB.human
    }
    else if(DatabaseType=='mouse'){
        CellChatDB<-CellChatDB.mouse
    }
    else{
        CellChatDB<-CellChatDB.human
    }

    showDatabaseCategory(CellChatDB)



    CellChatDB.use<-subsetDB(CellChatDB,search="Secreted Signaling")
    cellchat@DB<-CellChatDB.use

    cellchat<-subsetData(cellchat)

    cellchat <- identifyOverExpressedGenes(cellchat)

    cellchat <- identifyOverExpressedInteractions(cellchat)

    cellchat <- projectData(cellchat, PPI.human)


    cellchat <- computeCommunProb(cellchat)

    cellchat<-filterCommunication(cellchat,min.cells=10)
    df.net<-subsetCommunication(cellchat)

    # cellchat <- computeCommunProbPathway(cellchat) 
    # df.netp<-subsetCommunication(cellchat,slot.name="netP")



    ## 生成net
    cellchat<-aggregateNet(cellchat)

    ## 打包结果
    ### 通讯的count
    byCount <- list(cellchat@net$count,colnames(cellchat@net$count))
    names(byCount) <- c("data","clusterList")
    ### 通讯的weight
    byWeight <- list(cellchat@net$weight,colnames(cellchat@net$weight))
    names(byWeight) <- c("data","clusterList")
    ### Interaction的信息
    source <- df.net$source
    target <- df.net$target
    ligand <- df.net$ligand
    receptor <- df.net$receptor
    prob <- df.net$prob
    
    
    result <- list(byCount,byWeight,source,target,ligand,receptor,prob)
    names(result) <- c("count","weight","source","target","ligand","receptor","prob")

    return(result)

}


library(Matrix)
library(jsonlite)

# 保存的路径

args <- commandArgs(trailingOnly = TRUE)
JobId <- args[1]
save_path = file.path(".", "job_attachment", JobId, "CellChat")
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
result <- CellChat(matrix = mat, cellName = rownames, geneName = colnames, group = other_params$label, DatabaseType = other_params$organism)
write_json(result, file.path(save_path, "result.json"), pretty = TRUE, auto_unbox = TRUE)
