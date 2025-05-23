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
st = readRDS(file = paste0("../results/visium_axolotl_R17246_R12830_allVisium_20240905/Rdata/",
                           'seuratObject_allVisiusmst_',
                           'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                           '_R17246_R12830_allVisium_20240905', '.rds'))

Idents(st) = st$seg_ventricle
st = subset(st, cells = colnames(st)[which(!is.na(st$seg_ventricle))])

SpatialPlot(st, group.by = 'seg_ventricle', ncol = 4)

table(st$condition)

#InstallData("stxBrain")
#brain <- LoadData("stxBrain", type = "anterior1")
#load(file = paste0('../results/Rdata/', 
#                   'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))

st$condition = factor(st$condition)

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
  # n = 5
  cat(n, '  slice -- ', cc[n], '\n')
  slice = cc[n]
  
  outDir = paste0(resDir, '/', slice, '/')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  stx = subset(st, condition == slice)
  DefaultAssay(stx) = 'Spatial'
  
  # Normalization 
  stx <- SCTransform(stx, assay = "Spatial", verbose = FALSE)
  
  ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
  SpatialFeaturePlot(stx, features = ggs, images = cc[n], max.cutoff = 'q5')
  
  #ggsave(paste0(outDir, 'GAS6_AXL_expression_SCT_v2.pdf'), 
  #       width = 14, height = 6)
  
  # Dimensional reduction with all cells
  stx <- RunPCA(stx, assay = "SCT", verbose = FALSE)
  stx <- FindNeighbors(stx, reduction = "pca", dims = 1:30)
  stx <- FindClusters(stx, verbose = FALSE, resolution = 1.5)
  stx <- RunUMAP(stx, reduction = "pca", dims = 1:30)
  p1 <- DimPlot(stx, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
  p2 <- SpatialDimPlot(stx, label = TRUE, group.by = 'seurat_clusters', label.size = 3, images = cc[n])
  p1 + p2
  
  #ggsave(paste0(outDir, 'visium_clusters_v2.pdf'), width = 14, height = 6)
  
  
  ##########################################
  # # Format Spatial Coordinates and Normalize
  ##########################################
  stx@meta.data$x <- eval(parse(text = paste0("stx@images$", cc[n], "@coordinates$row")))
  stx@meta.data$y <- eval(parse(text = paste0("stx@images$", cc[n], "@coordinates$col")))
  
  DefaultAssay(stx) <- "Spatial"
  stx <- NormalizeData(stx)
  
  ggs = rownames(stx)[grep('GAS6|AXL', rownames(stx))]
  SpatialFeaturePlot(stx, features = ggs, images = cc[n])
  
  ggsave(paste0(outDir, 'GAS6_AXL_expression_logNormalize_v2.pdf'), 
         width = 14, height = 6)
  
  
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
    
    ggsave(paste0(outDir, 'GAS6_AXL_expression_ALRA.imputation.pdf'), 
           width = 14, height = 6)
    
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
                               #blend = 'sum',
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
    
    # # Find markers
    # mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T, test.use = "roc")
    # GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)
    # DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
    #   scale_fill_gradientn(colors = c("grey","white", "blue"))
    
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
    
    saveRDS(stx, file = paste0(outDir, 'st_res_NICHES_NeighborhoodToCell_alra.rds'))
    
    SpatialFeaturePlot(stx,
                       features = c('GAS6—AXL', "NRG1—ERBB2"),
                       slot = 'scale.data', images = slice)
    
    ggsave(paste0(outDir, '/GAS6_AXL_examples_interaction_visium_ALRAimputation.pdf'), 
           width = 14, height = 6)
    
  }
  
  ##########################################
  # use snRNA-seq data to imputate the visium data 
  ##########################################
  Use_snRNAseq_imputation = TRUE
  if(Use_snRNAseq_imputation){
    DefaultAssay(stx) <- "Spatial"
    # n = 3
    source('functions_Visium.R')
    RCTD_out = paste0("../results/visium_axolotl_R17246_R12830_allVisium_20240905/",
                      "RCTD_out/RCTD_subtype_out_41subtypes_ref.time.specific_v3.7_ventricleRegion/", 
                      slice)
    stx2 <- Run_imputation_snRNAseq_visium(stx, refs, 
                                           RCTD_out = RCTD_out,
                                           slice = slice, normalized_weights = FALSE)
    
    DefaultAssay(stx2) = 'imputated'
    ggs = rownames(stx2)[grep('GAS6|AXL|NRG1|ERBB2|AGRN|DAG1|POSTN|ITGB1', rownames(stx2))]
    SpatialFeaturePlot(stx2, features = ggs, images = slice)
    
    ggsave(paste0(outDir, 'GAS6_AXL_otherExamples_expression_visium_snRNAimputation.pdf'), 
           width = 14, height = 12)
    
    cat('change the feature names  \n')
    mat = stx2@assays$imputated@data
    meta = stx2@meta.data
    ggs = get_geneName(rownames(mat))
    mat = mat[match(unique(ggs), ggs), ]
    rownames(mat) = get_geneName(rownames(mat))
    
    ggs = rownames(mat)[grep('GAS6|AXL|NRG1|ERBB2|AGRN|DAG1|POSTN|ITGB1', rownames(mat))]
    
    srat <- CreateSeuratObject(counts = mat, data = mat,  assay = "imputation", 
                               meta.data = meta) # create object
    # 
    # image <- Read10X_Image(image.dir = paste0(topdir, "mock/outs/spatial/"),
    #                        filter.matrix = FALSE) # read in the images
    # image <- image[Cells(x = srat)] # filter image by the spots
    # DefaultAssay(object = image) <- "Spatial" # set default assay
    # srat[[keyname]] <- image # slice name might be changed
    # 
    NICHES_output <- RunNICHES(object = srat,
                               #LR.database = "fantom5",
                               LR.database = "custom",
                               species = "human",
                               assay = "imputation",
                               position.x = 'x',
                               position.y = 'y',
                               k = 4, 
                               #blend = 'sum',
                               custom_LR_database = data.frame(c('GAS6', 'NRG1', 'AGRN', 'POSTN'), 
                                                               c('AXL', 'ERBB2', 'DAG1', 'ITGB1')), 
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
    # mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T, test.use = "roc")
    # GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)
    # DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
    #   scale_fill_gradientn(colors = c("grey","white", "blue"))
    
    # Check that these make sense and print little plots
    # DefaultAssay(stx2) <- 'alra'
    # p1 <- SpatialFeaturePlot(stx2, crop = TRUE, features = "Fgf1",slot = "data",min.cutoff =  'q1',
    #                          max.cutoff = 'q99')+ggtitle("Ligand")+theme(legend.position = "right")
    # p2 <- SpatialFeaturePlot(stx2, crop = TRUE, features = "Fgfr2",slot = "data",min.cutoff =  'q1',
    #                          max.cutoff = 'q99')+ggtitle("Receptor")+theme(legend.position = "right")
    # 
    # ggpubr::ggarrange(p1,p2)
    
    # Add Niches output as an assay
    niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
    colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
    stx2[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
    
    
    #DefaultAssay(stx2) = 'alra'
    #ggs = rownames(stx2)[grep('GAS6|AXL', rownames(stx2))]
    #SpatialFeaturePlot(stx2, features = ggs, images = cc[n])
    
    # Plot celltype specific niche signaling
    DefaultAssay(stx2) <- "NeighborhoodToCell"
    stx2 <- ScaleData(stx2)
    saveRDS(stx2, file = paste0(outDir, 
                                'st_res_NICHES_NeighborhoodToCell_snRNA.imputation_LRexamples_batchAll.rds'))
    
    stx2 = readRDS(file = paste0(outDir, 'st_res_NICHES_NeighborhoodToCell_snRNA.imputation.rds'))
    
    SpatialFeaturePlot(stx2,
                       #features = c('GAS6—AXL', "NRG1—ERBB2", "AGRN-DAG1", "POSTN-ITGB1"),
                       features = rownames(stx2),
                       slot = 'scale.data', images = slice, 
                       min.cutoff = c(-1, -1),
                       max.cutoff = c(2, 3))
    
    ggsave(paste0(outDir, 'LRexamples_interaction_visium_snRNAimputation', cc[n], '_v2.pdf'), 
           width = 14, height = 12)
    
  
  }
  
}
