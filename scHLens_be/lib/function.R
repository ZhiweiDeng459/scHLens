
# A method of feature selection
# input:
    # matrix:non-negtive integer matrix(csc_matrix)
    # gene_name:a list of genes
    # cell_name:a list of cells
# output:a ranked list
devianceFS <- function(matrix,gene_name,cell_name){

    library(scry)

    rownames(matrix) = cell_name
    colnames(matrix) = gene_name

    rna_tmp <-devianceFeatureSelection(t(matrix))
    rna_tmp <- as.data.frame(rna_tmp)
    rna_tmp$gene_name <- rownames(rna_tmp)
    rna_tmp <- rna_tmp[order(-rna_tmp$rna_tmp),]
    scry_gene_list <- rna_tmp

    return(scry_gene_list$gene_name)

}

# A method of feature selection
# input:
    # matrix:non-negtive integer matrix(csc_matrix)
    # gene_name:a list of genes
    # cell_name:a list of cells
# output:a ranked list

MySCTransform <- function(matrix,gene_name,cell_name,topGenes){

    library(Seurat)

    rownames(matrix) = cell_name
    colnames(matrix) = gene_name

    rna <- CreateSeuratObject(counts = t(matrix))

    rna <- SCTransform(rna,variable.features.n=topGenes,verbose = FALSE,seed.use=1448145)
    
    result <- as.data.frame(rna@assays$SCT@scale.data)

    return(result)

}

# A method of trajectory inference
# input:
    # rd:the DR matrix
    # cl:a list of cell's annonation
# output: a S4 object

Slingshot <- function(rd,cl){

    library(SingleCellExperiment, quietly = TRUE)
    library(slingshot, quietly = TRUE)
    library(RColorBrewer)

    zeros = matrix(0,nrow(rd),2) ## fake count

    rownames(zeros) <- rownames(rd) ## prevent a bioconductor bug from unsame cell names

    sce <- SingleCellExperiment(assays=list(count=t(zeros)),reducedDims = SimpleList(rd=rd))



    colData(sce)$cl <- unlist(cl)


    ss_sce <- slingshot(sce, clusterLabels = 'cl', reducedDim = 'rd')

    

    slingshotList <- ss_sce@colData@listData



    ## generate colors
    PseudotimeValue <- slingshotList[grep('^slingPseudotime_.*',names(slingshotList))]
    PseudotimeValue <- data.frame(PseudotimeValue)

    PseudotimeAgg <- apply(PseudotimeValue,1,max,na.rm=T)


    colors <- colorRampPalette(brewer.pal(11,'Spectral'))(100)

    PseudotimeColor <- colors[cut(PseudotimeAgg, breaks=100)]
                    
    ## generate curves
    lin <- slingshotList$slingshot@metadata$lineages
    crv <- slingshotList$slingshot@metadata$curves

    ## pack
    result <- list(PseudotimeColor=PseudotimeColor,lineages = lin,curves = crv)

}


# A method of cell chat
# input:
    # matrix: the normlized matrix
    # cellName: A list of cell name
    # geneName: A list of gene name
    # group: a list of cell's annonation
    # DatabaseType: human or mouse
# output: an ListVector

CellChat <- function(matrix,cellName,geneName,group,DatabaseType){

    library(scater)
    library(Seurat)
    library(SeuratDisk)
    library(patchwork)
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

    cellchat <- computeCommunProbPathway(cellchat) 
    df.netp<-subsetCommunication(cellchat,slot.name="netP")



    ## 生成net
    cellchat<-aggregateNet(cellchat)

    ## 打包结果
    byCount <- list(cellchat@net$count,colnames(cellchat@net$count))
    names(byCount) <- c("data","clusterList")
    byWeight <- list(cellchat@net$weight,colnames(cellchat@net$weight))
    names(byWeight) <- c("data","clusterList")
    result <- list(byCount,byWeight)
    names(result) <- c("count","weight") 

    return(result)

}