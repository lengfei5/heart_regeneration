##########################################################################
##########################################################################
# Project: heart regeneration project
# Script purpose: First script to analyze the processed Visium data by spaceranger of 10x
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Aug 27 11:39:22 2021
##########################################################################
##########################################################################
# setup for data import and sequencing QCs
version.analysis = '_R11934_20210827'

resDir = paste0("../results/visium_mouse", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium'

require(Seurat)
require(SeuratObject)
require(ggplot2)
require(tibble)
require(dplyr)
library(patchwork)

########################################################
########################################################
# Section : import the processed visium data by spaceranger
# first start with neonatal mice samples
########################################################
########################################################
#design = data.frame(seq(166904, 166911), c(paste0("adult.day", c(14, 7, 4, 1)), 
#                                      paste0('neonatal.day', c(1, 4, 7, 14))), stringsAsFactors = FALSE)
design = data.frame(seq(166908, 166911), c(paste0('neonatal.day', c(1, 4, 7, 14))), stringsAsFactors = FALSE)
design = data.frame(seq(166904, 166907), c(paste0("adult.day", c(14, 7, 4, 1))), stringsAsFactors = FALSE)
colnames(design) = c('sampleID', 'condition')

varibleGenes = c()
for(n in 1:nrow(design))
{
  # n = 4
  
  # load output from spaceranger
  aa = Seurat::Load10X_Spatial(
    data.dir = paste0(dataDir, '/output_', design$sampleID[n], '/', design$sampleID[n],  '/outs'),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice =  design$condition[n],
    filter.matrix = TRUE,
    to.upper = FALSE
  )
  
  aa$condition = design$condition[n]
  #aa <- SCTransform(aa, assay = "Spatial",  method = "glmGamPoi", verbose = FALSE)
  aa <- SCTransform(aa, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE)
  
  test.clustering.each.condtiion = FALSE
  if(test.clustering.each.condtiion){
    aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(aa)
    
    aa <- FindNeighbors(aa, dims = 1:10)
    aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
    aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.05)
    
    DimPlot(aa, reduction = "umap", group.by = c("ident"))
    
  }
  
  varibleGenes = unique(c(varibleGenes, VariableFeatures(aa)))
  
  cat(design$condition[n], ' : ',  ncol(aa), ' spot found \n')
  
  # merge slices from different time points and 
  if(n == 1) {
    st = aa
  }else{
    st = merge(st, aa)
  }
  
}


##########################################
# cell and gene filtering
##########################################
#st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^Mt-")
# Visualize QC metrics as a violin plot
VlnPlot(st, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 2)

Idents(st) = st$condition
FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

#st = SCTransform(st, assay = "Spatial", verbose = FALSE)

DefaultAssay(st) <- "SCT"

VariableFeatures(st) <- varibleGenes

st <- RunPCA(st, verbose = FALSE)
st <- FindNeighbors(st, dims = 1:30)
st <- FindClusters(st, verbose = FALSE)
st <- RunUMAP(st, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(st, reduction = "umap", group.by = c("ident", "condition"))


#plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
#plot2 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "right")
#wrap_plots(plot1, plot2)

SpatialDimPlot(st)

SpatialFeaturePlot(st, features = c("Myh6"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Tcf21"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Vim"), image.alpha = 0.5)

SpatialFeaturePlot(st, features = c("Ddr2"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Acta2"), image.alpha = 0.5)

SpatialFeaturePlot(st, features = c("Emcn"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Kdr"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Ptprc"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Cd68"), image.alpha = 0.5)

SpatialFeaturePlot(st, features = c("Itgam"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c(""), image.alpha = 0.5)

SpatialFeaturePlot(st, features = c("Fstl1"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Wt1"), image.alpha = 0.5)
SpatialFeaturePlot(st, features = c("Agrn"), image.alpha = 0.5)
