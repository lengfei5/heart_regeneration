rm(list = ls())

version.analysis = '_scRNAseq_20250227'

resDir = paste0("../results/scRNAseq_zebrafish", version.analysis,'/')
RdataDir = paste0(resDir, 'Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

# required libraries
library(data.table)
require(Seurat)
library(SeuratObject)
require(sctransform)
require(ggplot2)
library(dplyr)
library(patchwork)
require(tictoc)
library(pryr) # monitor the memory usage
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(reshape2)


source('functions_scRNAseq.R')
source('functions_Visium.R')
source('utility_zebrafish.R')

options(future.globals.maxSize = 32000 * 1024^2)

mem_used()

dataDir = "../published_dataset/zebrafish/Hu_Junker_2022/zebrafish_heart_processing/"


source('functions_dataIntegration.R')

outDir = paste0(resDir, '/batch_correction/')
if(!dir.exists(outDir)) dir.create(outDir)

aa = readRDS(file = paste0(RdataDir, 'dr_scRNAseq_39Batches_noCorrection.rds'))

aa  <-  aa %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = 'pca', n.neighbors = 30, min.dist = 0.1)

## test Harmony
aa <- harmony::RunHarmony(aa, 
                          group.by.vars = c("batch"),
                          reduction = "pca", 
                          max.iter.harmony = 30, 
                          dims.use = 1:30,
                          project.dim = FALSE,
                          plot_convergence = TRUE)

aa <- RunUMAP(aa, dims = 1:30, 
                reduction = 'harmony',
                reduction.name = "umap_harmony",
                reduction.ke = 'umapharmony_',
                verbose = FALSE,
                min.dist = 0.1)

saveRDS(aa, file = paste0(RdataDir, '/dr_scRNAseq_39Batches_correction_Harmony.rds'))

# p1 = DimPlot(aa, group.by = "celltypes", label = TRUE, repel = TRUE, raster=FALSE, 
#              reduction = 'umap_harmony')
# p2 = DimPlot(aa, group.by = "batch", label = FALSE, raster=FALSE, reduction = 'umap_harmony')
# 
# p1 / p2

ggsave(filename = paste0(resDir, 'UMAP_celltypes.original_39batch_dataIntegration_Harmony.pdf'), 
       height = 16, width = 12)


source('functions_dataIntegration.R')
ref.combined = IntegrateData_runFastMNN(aa, group.by = 'batch', 
                                        nfeatures = 3000,
                                        #merge.order = list(list(4,3,1,2), list(8,7,5,6)),
                                        correct.all = TRUE)
DefaultAssay(ref.combined) = 'mnn.reconstructed'

saveRDS(ref.combined, file = paste0(outDir, 
                                    'dr_scRNAseq_39Batches_correction_fastMNN.rds'))


ref.combined<- RunUMAP(ref.combined, dims = 1:30, 
              reduction = 'mnn',
              reduction.name = "umap_mnn",
              reduction.ke = 'umapmnn_',
              verbose = FALSE,
              min.dist = 0.1)

p1 = DimPlot(aa, group.by = "celltypes", label = TRUE, repel = TRUE, raster=FALSE,
             reduction = 'umap_mnn')
p2 = DimPlot(aa, group.by = "batch", label = FALSE, raster=FALSE, reduction = 'umap_mnn')

p1 / p2
ggsave(filename = paste0(resDir, 'UMAP_celltypes.original_39batch_dataIntegration_mnn.pdf'), 
       height = 16, width = 12)



## RPCA integration
method = "Seurat_RPCA"
ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                         group.by = 'batch', 
                                         nfeatures = 3000,
                                         #merge.order = 
                                         redo.normalization.scaling = TRUE,
                                         correct.all = TRUE, 
                                         use.parallelization = TRUE)

DefaultAssay(ref.combined) = 'integrated'
saveRDS(ref.combined, file = paste0(outDir, 
                                    'dr_scRNAseq_39Batches_correction_SeuratRPCA.rds'))


