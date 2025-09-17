install.packages("devtools",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 2.4.5
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 1.30.25
BiocManager::install("S4Vectors") # 0.44.0
BiocManager::install("SingleCellExperiment") # 1.28.1
BiocManager::install("Seurat") # 5.2.1
BiocManager::install("SeuratDisk") # -
BiocManager::install("scater") # -
BiocManager::install("Signac") # 1.14.0
BiocManager::install("scry") # 1.18.0
BiocManager::install("dplyr") # 1.1.4
BiocManager::install("patchwork") # 1.3.0
BiocManager::install("reticulate") # 1.41.0.1
BiocManager::install("glmGamPoi") # 1.18.0
BiocManager::install("kstreet13/slingshot") # -
BiocManager::install("BiocNeighbors") # 2.0.1
BiocManager::install("ComplexHeatmap") # 2.22.0
BiocManager::install("limma") # 3.62.2
devtools::install_github('satijalab/seurat-data') # -
devtools::install_github("sqjin/CellChat") # -
devtools::install_github("saeyslab/nichenetr") # -
install.packages('Matrix',repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 1.7.3
install.packages("uwot",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 0.2.3
install.packages("mclust",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 6.1.1
install.packages("CellChat",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # -
install.packages("tidyverse",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 2.0.0
remotes::install_github("mojaveazure/seurat-disk")
install.packages(c("tensorflow", "keras"),repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE)
devtools::install_github('PYangLab/scCCESS')
devtools::install_github("yulijia/SIMLR", ref = 'master')
