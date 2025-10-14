#' MultiK main algorithm
#'
#' MultiK main algorithm: takes a preprocessed gene expression matrix as input. Then subsamples 80\% of the cells and applies standard Seurat pipeline on the subsampled data matrix 100 times over multiple resolution parameters.
#' @param seu A Seurat object with normalized count
#' @param resolution A vector Seurat resolution parameters. Default is from 0.05 to 2 with step size of 0.05
#' @param nPC Number of principal components to use in clustering
#' @param reps Number of subsampling runs. Integer value. Default is 100
#' @param pSample Proportion of cells to sample. Numerical value. Default is 0.8
#' @param seed Optional numerical value. This sets a random seed for generating reproducible results
#' @return A list with components: k is a vector of number of runs for each K. clusters is a list containing the clustering labels for each subsampling run at each resolution parameter. consensus is a list containing a consensus matrix for each K.
#' @export
MultiK1 <- function(seu, resolution = seq(0.05, 2, 0.05), nPC = 30, reps = 100, pSample = 0.8, seed = NULL) {
  # setting seed for reproducibility
  if (is.null(seed) == TRUE) {
    seed <- timeSeed <- as.numeric(Sys.time())
  }
  set.seed(seed)

  # step 1: subsampling
  subcol <- list()
  for (i in 1: reps) {
    subcol[[i]] <- sample(x=ncol(seu), size=round(ncol(seu) * pSample), replace=FALSE)
  }

  # step 2: loop over subsampling runs, with each run subsampling 80% of cells, reselect genes for clustering
  clusters <- list()
  messages <- c()
  ks <- c()
  count <- 1

  suppressPackageStartupMessages(library(Seurat))

  for(i in 1: reps) {

    print(paste("Rep: ", i, sep=""))
    # subsample the columns (the cells) from the full matrix
    subX <- seu[, subcol[[i]] ]

    # normalizing the data
    #subX <- NormalizeData(object = subX, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

    # Find HVG genes ~ 2000 genes
    subX <- FindVariableFeatures(object = subX, selection.method = "vst", nfeatures = 2000,
                                 loess.span = 0.3, clip.max = "auto",
                                 num.bin = 20, binning.method = "equal_width", verbose = F)

    # Scaling unwanted variation
    all.genes <- rownames(x = subX)
    subX <- ScaleData(object = subX, features = all.genes, verbose=F)
    # Run PCA to reduce dimensions
    subX <- RunPCA(object = subX, features = VariableFeatures(object = subX), npcs = 50, verbose=F)
    # Run Clustering
    subX <- FindNeighbors(object = subX,
                          k.param = 20, # default is 20-nearest neighbors
                          reduction = "pca", dims = 1: nPC, verbose=F)

    for (res in resolution) {
      print(paste("Rep", i, "Res", res, sep=" "))
      subX <- FindClusters(subX, resolution = res, verbose = F)
      subX.clusters <- Idents(subX)

      # skip length(unique(subX.clusters)) < 3 or > 15
      num_clusters <- length(unique(subX.clusters))
      if (num_clusters < 3 | num_clusters > 15) {
        print(paste("Skipping Rep", i, "Res", res, "due to cluster count:", num_clusters, sep=" "))
        next
      }

      clusters[[count]] <- subX.clusters
      messages <- c(messages, paste("Rep_", i, "_res_", res, sep = ""))
      count <- count + 1
      ks <- c(ks, length(unique(subX.clusters)))
    }
    names(clusters) <- messages

  }

  # step 3: calculate consensus matrix across subsampling runs for each unique K
  mInit <- matrix(0, ncol = ncol(seu), nrow = ncol(seu))

  ml <- list()
  res <- list()
  all.clusters.by.K <- list()
  m.count <- list()
  unique.ks <- unique(ks)[order(unique(ks))]

  count.k <- 1
  for(k in unique.ks) {
    print(paste("k =", k, sep=" "))
    idx <- which(ks == k)
    cluster.k <- clusters[idx]
    all.clusters.by.K[[count.k]] <- cluster.k

    for (s in 1: length(cluster.k) ) {
      print(paste("run", s, sep = ""))
      sampleKey <- as.numeric(sapply(names(cluster.k[[s]]), function(x){which(colnames(seu) == x)}))
      if (s == 1){
        ml[[count.k]] <- connectivityMatrix(cluster.k[[s]], mInit, sampleKey)
        m.count[[count.k]] <- connectivityMatrix(rep(1, length(sampleKey)), mInit, sampleKey)
      }else{
        ml[[count.k]] <- connectivityMatrix(cluster.k[[s]], ml[[count.k]], sampleKey)
        m.count[[count.k]] <- connectivityMatrix(rep(1, length(sampleKey)), m.count[[count.k]], sampleKey)
      }
    }

    res[[count.k]] <- triangle(ml[[count.k]], mode = 3)/triangle(m.count[[count.k]], mode = 3)
    res[[count.k]][which(triangle(m.count[[count.k]], mode = 3) == 0)] = 0
    print(paste(k, " finished", sep = ""))
    count.k <- count.k + 1
  }

  return(list("consensus" = res, "k" = ks))
}


multik <- function(norm_data,gene_name,cell_name){

    # norm_data: normalizad matrix
    # gene_name: a gene name list
    # cell_name: a cell name list

    library(MultiK)
    library(Seurat)

    rownames(norm_data) = cell_name
    colnames(norm_data) = gene_name
    norm_data <- t(norm_data)

    seu <- CreateSeuratObject(counts=norm_data) 
    seu <- SetAssayData(seu, layer = "data", new.data = norm_data)


    rep_num = 10

    multik <- MultiK1(seu, resolution = seq(0.01, 1, 0.05), reps=rep_num,seed = 42)

    ## calc 1-rPCA (x) / runs (y)
    tog <- as.data.frame(table(multik$k)[table(multik$k) > 1])
    pacobj <- CalcPAC(x1=0.1, x2=0.9, xvec = tog$Var1, ml = multik$consensus)
    tog$rpac <- pacobj$rPAC
    tog$one_minus_rpac  <- 1-tog$rpac


    hull_idx <- chull(tog$one_minus_rpac,tog$Freq)

    raw_hull_points <- tog[hull_idx, ]


    hull_points <- raw_hull_points[raw_hull_points$Freq > (1 / nrow(tog)) * rep_num * 20,]
    if(nrow(hull_points)==0){
        hull_points <- raw_hull_points
    }

    best_k = hull_points[which.max(hull_points$one_minus_rpac),'Var1']

    return(best_k)

}

library(Matrix)
library(jsonlite)

# 保存的路径
args <- commandArgs(trailingOnly = TRUE)
JobId <- args[1]
save_path = file.path(".", "job_attachment", JobId, "multik")

# 读取矩阵
mat <- readMM(file.path(save_path, "matrix.mtx"))
mat <- as(mat, "CsparseMatrix")
  
# 读取行名、列名
rownames <- read.table(file.path(save_path, "rows.tsv"), stringsAsFactors = FALSE)[,1]
colnames <- read.table(file.path(save_path, "cols.tsv"), stringsAsFactors = FALSE)[,1]


# 如果你的 multik 函数定义在 R 里，可以直接调用
result <- multik(norm_data = mat, gene_name = colnames, cell_name = rownames)

result_json <- list(
  best_k = result
)
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
write_json(result_json, file.path(save_path, "result.json"), pretty = TRUE)
