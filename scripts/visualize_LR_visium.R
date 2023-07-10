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
source('functions_scRNAseq.R')
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

refs = readRDS(file = paste0(RdataDir, 'RCTD_refs_subtypes_final_20221117.rds'))
refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$celltypes = refs$celltype_toUse

table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = as.factor(refs$subtypes) 
refs$celltypes = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltypes)
refs$celltypes = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', refs$celltypes)


##########################################
# we need to loop over the time points
##########################################
cat('-- check visium conditions -- \n')
st$condition = droplevels(as.factor(st$condition))
print(table(st$condition))
cc = names(table(st$condition))

Idents(st) = st$condition

for(n in 1:length(cc))
{
  # n = 3
  cat(n, '  slice -- ', cc[n], '\n')
  slice = cc[n]
  stx = subset(st, condition == slice)
  DefaultAssay(stx) = 'Spatial'
  
  # Normalization 
  stx <- SCTransform(stx, assay = "Spatial", verbose = FALSE)
  
  ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
  SpatialFeaturePlot(stx, features = ggs, images = cc[n])
  
  # Dimensional reduction with all cells
  stx <- RunPCA(stx, assay = "SCT", verbose = FALSE)
  stx <- FindNeighbors(stx, reduction = "pca", dims = 1:30)
  stx <- FindClusters(stx, verbose = FALSE, resolution = 1.0)
  stx <- RunUMAP(stx, reduction = "pca", dims = 1:30)
  p1 <- DimPlot(stx, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
  p2 <- SpatialDimPlot(stx, label = TRUE, group.by = 'seurat_clusters', label.size = 3, images = cc[n])
  p1 + p2

  ##########################################
  # # Format Spatial Coordinates and Normalize
  ##########################################
  stx@meta.data$x <- eval(parse(text = paste0("stx@images$", cc[n], "@coordinates$row")))
  stx@meta.data$y <- eval(parse(text = paste0("stx@images$", cc[n], "@coordinates$col")))
  
  DefaultAssay(stx) <- "Spatial"
  stx <- NormalizeData(stx)
  
  ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
  SpatialFeaturePlot(stx, features = ggs, images = cc[n])
  
  ##########################################
  # # Impute and Run NICHES
  ##########################################
  Use_ALRA_impuation = FALSE
  if(Use_ALRA_impuation){
    DefaultAssay(stx) <- "Spatial"
    stx <- SeuratWrappers::RunALRA(stx)
    
    DefaultAssay(stx) = 'alra'
    ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
    SpatialFeaturePlot(stx, features = ggs, images = cc[n])
    
    cat('change the feature names  \n')
    mat = stx@assays$alra@data
    meta = stx@meta.data
    ggs = get_geneName(rownames(mat))
    mat = mat[match(unique(ggs), ggs), ]
    rownames(mat) = get_geneName(rownames(mat))
    
    srat <- CreateSeuratObject(counts = mat, data = mat,  assay = "alra", meta.data = meta) # create object
    # 
    # image <- Read10X_Image(image.dir = paste0(topdir, "mock/outs/spatial/"),
    #                        filter.matrix = FALSE) # read in the images
    # image <- image[Cells(x = srat)] # filter image by the spots
    # DefaultAssay(object = image) <- "Spatial" # set default assay
    # srat[[keyname]] <- image # slice name might be changed
    # 
    NICHES_output <- RunNICHES(object = srat,
                               LR.database = "fantom5",
                               species = "human",
                               assay = "alra",
                               position.x = 'x',
                               position.y = 'y',
                               k = 4, 
                               blend = 'sum',
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
    # DefaultAssay(stx) <- 'alra'
    # p1 <- SpatialFeaturePlot(stx, crop = TRUE, features = "Fgf1",slot = "data",min.cutoff =  'q1',
    #                          max.cutoff = 'q99')+ggtitle("Ligand")+theme(legend.position = "right")
    # p2 <- SpatialFeaturePlot(stx, crop = TRUE, features = "Fgfr2",slot = "data",min.cutoff =  'q1',
    #                          max.cutoff = 'q99')+ggtitle("Receptor")+theme(legend.position = "right")
    # 
    # ggpubr::ggarrange(p1,p2)
    
    # Add Niches output as an assay
    niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
    colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
    stx[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
    
    
    DefaultAssay(stx) = 'alra'
    ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
    SpatialFeaturePlot(stx, features = ggs, images = cc[n])
    
    ggsave(paste0(resDir, '/GAS6_AXL_expression_visium_ALRAimputation', cc[n], '.pdf'), 
           width = 10, height = 6)
    
    
    # Plot celltype specific niche signaling
    DefaultAssay(stx) <- "NeighborhoodToCell"
    stx <- ScaleData(stx)
    SpatialFeaturePlot(stx,
                       features = c('GAS6—AXL'),
                       slot = 'scale.data', images = cc[n])
    
    ggsave(paste0(resDir, '/GAS6_AXL_interaction_visium_ALRAimputation', cc[n], '.pdf'), 
           width = 10, height = 6)
    
    
  }
  
  ##########################################
  # use snRNA-seq data to imputate the visium data 
  ##########################################
  Use_snRNAseq_imputation = FALSE
  if(Use_snRNAseq_imputation){
    DefaultAssay(stx) <- "Spatial"
    # n = 3
    source('functions_Visium.R')
    slice = cc[n]
    stx <- Run_imputation_snRNAseq_visium(stx, refs, slice = slice, normalized_weights = FALSE)
    
    DefaultAssay(stx) = 'imputated'
    ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
    SpatialFeaturePlot(stx, features = ggs, images = slice)
    
    ggsave(paste0(resDir, '/GAS6_AXL_expression_visium_snRNAimputation', cc[n], '.pdf'), 
           width = 10, height = 6)
    
    cat('change the feature names  \n')
    mat = stx@assays$imputated@data
    meta = stx@meta.data
    ggs = get_geneName(rownames(mat))
    mat = mat[match(unique(ggs), ggs), ]
    rownames(mat) = get_geneName(rownames(mat))
    
    srat <- CreateSeuratObject(counts = mat, data = mat,  assay = "imputation", meta.data = meta) # create object
    # 
    # image <- Read10X_Image(image.dir = paste0(topdir, "mock/outs/spatial/"),
    #                        filter.matrix = FALSE) # read in the images
    # image <- image[Cells(x = srat)] # filter image by the spots
    # DefaultAssay(object = image) <- "Spatial" # set default assay
    # srat[[keyname]] <- image # slice name might be changed
    # 
    NICHES_output <- RunNICHES(object = srat,
                               LR.database = "fantom5",
                               species = "human",
                               assay = "imputation",
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
    # DefaultAssay(stx) <- 'alra'
    # p1 <- SpatialFeaturePlot(stx, crop = TRUE, features = "Fgf1",slot = "data",min.cutoff =  'q1',
    #                          max.cutoff = 'q99')+ggtitle("Ligand")+theme(legend.position = "right")
    # p2 <- SpatialFeaturePlot(stx, crop = TRUE, features = "Fgfr2",slot = "data",min.cutoff =  'q1',
    #                          max.cutoff = 'q99')+ggtitle("Receptor")+theme(legend.position = "right")
    # 
    # ggpubr::ggarrange(p1,p2)
    
    # Add Niches output as an assay
    niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
    colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
    stx[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
    
    
    #DefaultAssay(stx) = 'alra'
    #ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
    #SpatialFeaturePlot(stx, features = ggs, images = cc[n])
    
    # Plot celltype specific niche signaling
    DefaultAssay(stx) <- "NeighborhoodToCell"
    stx <- ScaleData(stx)
    SpatialFeaturePlot(stx,
                       features = c('GAS6—AXL'),
                       slot = 'scale.data', images = cc[n])
    
    ggsave(paste0(resDir, '/GAS6_AXL_interaction_visium_snRNAimputation', cc[n], '.pdf'), 
           width = 10, height = 6)
    
  }
  
  
}
