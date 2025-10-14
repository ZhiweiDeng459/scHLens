# install.packages("devtools",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 2.4.5
install.packages("BiocManager",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 1.30.25
install.packages("remotes",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE)
remotes::install_version("BiocManager", version = "1.30.19")
# options("repos"=c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

BiocManager::install("S4Vectors") # 0.44.0
BiocManager::install("SingleCellExperiment") # 1.28.1
remotes::install_version("Seurat", version = "5.2") # 5.2.1
# BiocManager::install("SeuratDisk") # 无
BiocManager::install("scater") # 无
BiocManager::install("Signac") # 1.14.0
BiocManager::install("GenomeInfoDbData")
BiocManager::install("scry") # 1.18.0
BiocManager::install("dplyr") # 1.1.4
BiocManager::install("patchwork") # 1.3.0
BiocManager::install("reticulate") # 1.41.0.1
BiocManager::install("glmGamPoi") # 1.18.0
# BiocManager::install("kstreet13/slingshot") # 无
BiocManager::install("BiocNeighbors") # 2.0.1
BiocManager::install("ComplexHeatmap") # 2.22.0
BiocManager::install("limma") # 3.62.2
# remotes::install_github('satijalab/seurat-data') # 无
remotes::install_github("sqjin/CellChat") # 
# remotes::install_github("saeyslab/nichenetr")
# install.packages("uwot",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 0.2.3
install.packages("mclust",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 6.1.1
# install.packages("CellChat",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 无
install.packages("tidyverse",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # 2.0.0
remotes::install_github("siyao-liu/MultiK")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages(c("tensorflow", "keras"),repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE)
# remotes::install_github('PYangLab/scCCESS')
# remotes::install_github("yulijia/SIMLR", ref = 'master')
# install.packages("plumber",repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",quietly=TRUE) # TODO

