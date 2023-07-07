##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: visualize LR interaction in the visium data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jul  7 12:22:45 2023
##########################################################################
##########################################################################

########################################################
########################################################
# Section : test NICHES
# #  https://msraredon.github.io/NICHES/articles/01%20NICHES%20Spatial.html 
########################################################
########################################################
rm(list = ls())

library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(NICHES)
library(viridis)


species = 'axolotl'
version.analysis = '_R12830_resequenced_20220308'
resDir = paste0("../results/visium_axolotl", version.analysis, '/LR_visualization')
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
options(future.globals.maxSize = 100000 * 1024^2)

mem_used()

##########################################
# Load Data, Normalize, Visualize
##########################################
#InstallData("stxBrain")
#brain <- LoadData("stxBrain", type = "anterior1")
load(file = paste0('../results/Rdata/', 
                   'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
st$condition = factor(st$condition, levels = design$condition)


# Normalization 
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

# Dimensional reduction with all cells
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE,group.by = 'seurat_clusters', label.size = 3)
p1 + p2

##########################################
# # Format Spatial Coordinates and Normalize
##########################################
brain@meta.data$x <- brain@images$anterior1@coordinates$row
brain@meta.data$y <- brain@images$anterior1@coordinates$col

DefaultAssay(brain) <- "Spatial"
brain <- NormalizeData(brain)

##########################################
# # Impute and Run NICHES
##########################################
brain <- SeuratWrappers::RunALRA(brain)

NICHES_output <- RunNICHES(object = brain,
                           LR.database = "fantom5",
                           species = "mouse",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           k = 4, 
                           cell_types = "seurat_clusters",
                           min.cells.per.ident = 0,
                           min.cells.per.gene = NULL,
                           meta.data.to.map = c('orig.ident','seurat_clusters'),
                           CellToCell = F, CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F, 
                           NeighborhoodToCell = T)


niche <- NICHES_output[['NeighborhoodToCell']]
Idents(niche) <- niche[['ReceivingType']]

# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)

niche <- RunUMAP(niche,dims = 1:10)
DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +
  ggtitle('Cellular Microenvironment') + 
  NoLegend()

# Find markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T, test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)
DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
  scale_fill_gradientn(colors = c("grey","white", "blue"))

# Check that these make sense and print little plots
DefaultAssay(brain) <- 'alra'
p1 <- SpatialFeaturePlot(brain, crop = TRUE, features = "Fgf1",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99')+ggtitle("Ligand")+theme(legend.position = "right")
p2 <- SpatialFeaturePlot(brain, crop = TRUE, features = "Fgfr2",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99')+ggtitle("Receptor")+theme(legend.position = "right")

ggpubr::ggarrange(p1,p2)


# Add Niches output as an assay
niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
brain[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
DefaultAssay(brain) <- "NeighborhoodToCell"
brain <- ScaleData(brain)


# Plot celltype specific niche signaling
SpatialFeaturePlot(brain,
                   features = c('Bmp2—Bmpr2','Efna1—Ephb6','Fgf1—Fgfr2'),
                   slot = 'scale.data')



