rm(list = ls())

version.analysis = '_R13591_intron.exon.20220729'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome'

source('functions_scRNAseq.R')
source('functions_Visium.R')

require(Seurat)
#require(sctransform)
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 300000 * 1024^2)

mem_used()

species = 'axloltl_scRNAseq'

##########################################
# main code
##########################################
sub.obj = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                                '_lognormamlized_pca_umap_keep.missed_subtypes_DFinderFiltered_celltypes.rds'))
#sub.obj = aa
#rm(aa)

sub.list <- SplitObject(sub.obj, split.by = "condition")

# normalize and identify variable features for each dataset independently
sub.list <- lapply(X = sub.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = sub.list, nfeatures = 3000)
features.common = rownames(sub.obj)
sub.list <- lapply(X = sub.list, FUN = function(x) {
  x <- ScaleData(x, features = features.common, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  
})

sub.anchors <- FindIntegrationAnchors(object.list = sub.list, 
                                      anchor.features = features, 
                                      reduction = "rpca", 
                                      k.anchor = 5)

rm(sub.list)

# this command creates an 'integrated' data assay
sub.combined <- IntegrateData(anchorset = sub.anchors, features.to.integrate = features.common)

rm(sub.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(sub.combined) <- "integrated"

xx = DietSeurat(sub.combined, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'integrated')
xx@assays$integrated@counts = sub.combined@assays$RNA@counts


# Run the standard workflow for visualization and clustering
sub.combined = xx;
rm(xx)
sub.combined$condition = factor(sub.combined$condition, 
                                levels = c('Amex_scRNA_d0', 'Amex_scRNA_d1', 
                                           'Amex_scRNA_d4', 'Amex_scRNA_d7',
                                           'Amex_scRNA_d14'))

sub.combined <- ScaleData(sub.combined, verbose = FALSE)
sub.combined <- RunPCA(sub.combined, npcs = 30, verbose = FALSE)

#ElbowPlot(sub.combined, ndims = 30)

sub.combined <- FindNeighbors(sub.combined, reduction = "pca", dims = 1:30)
sub.combined <- FindClusters(sub.combined, algorithm = 3, resolution = 0.5)

sub.combined <- RunUMAP(sub.combined, reduction = "pca", dims = 1:30, n.neighbors = 30, min.dist = 0.3) 

#p0 = DimPlot(sub.obj, reduction = "umap")
#p1 = DimPlot(sub.combined, reduction = "umap", group.by = 'RNA_snn_res.0.5')
#p0 | p1

#DimPlot(sub.combined, reduction = 'umap', label = TRUE, label.size = 4,  split.by = 'condition') 
saveRDS(sub.combined, file = paste0(RdataDir, 'seuratObject_', species, 
                                   version.analysis, 
                                   '_lognormamlized_pca_umap_keep.missed_subtypes_DFinderFiltered_',
                                   'celltypes_timeBC_Rscript.rds'))
