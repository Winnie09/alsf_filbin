<<<<<<< HEAD
setwd('/users/whou/')
=======
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
>>>>>>> 493eb965c5992a7a722bc79e959aa16119876f2b
library(Seurat)
af = list.files('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/')
i = as.numeric(commandArgs(trailingOnly = T)[[1]])
d.integrated = readRDS(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/integrated.rds'))
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(d.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
d.integrated <- ScaleData(d.integrated, verbose = FALSE)
d.integrated <- RunPCA(d.integrated, npcs = 30, verbose = FALSE)
d.integrated <- RunUMAP(d.integrated, reduction = "pca", dims = 1:30)
u = d.integrated@reductions$umap@cell.embeddings
dir.create(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/umap/'),showWarnings = F,recursive = T)
saveRDS(u,paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/umap/umap.rds'))

pca = d.integrated@reductions$pca@cell.embeddings
dir.create(paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/pca/'),showWarnings = F,recursive = T)
saveRDS(pca,paste0('./alsf_filbin/data/integrate_each_to_atlas_rmDEG/result/',af[i],'/pca/pca.rds'))
