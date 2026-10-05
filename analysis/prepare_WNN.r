# To integrate ATAC-seq data and then perform multi-modal analyess with WNN

library(Seurat)
library(Signac)
#library(hdf5r) #need to load `hdf5_hl` on rcc
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(data.table)
library(ggplot2)

#immune.combined <- readRDS("~/cluster/projects/u19_multiomics/analyses/remove_dup/lungs_spleens_integrated_sub.RDS")

#DefaultAssay(immune.combined) <- "ATAC"
# atac.list <- SplitObject(immune.combined, split.by = "full.ident")
# atac.list <- lapply(X = atac.list, FUN = function(x) {
#     x <- RunTFIDF(x)
#     x <- FindTopFeatures(x, min.cutoff = 'q0')
#     x <- RunSVD(x)
#     x <- RunUMAP(x, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# })

atac.list <- readRDS("~/cluster/projects/u19_multiomics/analyses/remove_dup/atac_list_sub.RDS")

anchors <- FindIntegrationAnchors(object.list = atac.list, reduction = "rpca", dims = 1:20)
atac.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
atac.integrated <- ScaleData(atac.integrated, verbose = FALSE)
atac.integrated <- RunPCA(atac.integrated, verbose = FALSE)
atac.integrated <- RunUMAP(atac.integrated, dims = 1:20)
atac.integrated <- FindNeighbors(atac.integrated, reduction = "pca", dims = 1:20)
atac.integrated <- FindClusters(atac.integrated, resolution = 0.5)

saveRDS(atac.integrated, "~/cluster/projects/u19_multiomics/analyses/remove_dup/atac_integrated.RDS")