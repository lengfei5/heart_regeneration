##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: process the scRNA-seq data and ST data from Li et al. 2025
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jul  2 12:26:49 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_scRNAseq_Stereoseq_20260702'

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
mem_used()

dataDir = "../published_dataset/zebrafish/Li_et_al_2025/"


#mito.genes names
#mito.genes <- read.table("Data/mito.genes.vs.txt",sep = ",")
#mito.genes <- mito.genes$V3
#mito.genes <- as.character(mito.genes)

########################################################
########################################################
# Section I: import the data and prepare the scRNA-seq data
# 
########################################################
########################################################
## load and integrate data in one seurat object, calculate mito reads, filter cells, 

## the chamber scRNAseq data is from the Fig 5 in the original paper to construct a 3D map
#aa = readRDS(file = paste0(dataDir, 'scRNA-seq-chamber.rds'))

## the 206719 cells of regeneration is from Fig 1, the main dataset of regeneration
aa = readRDS(file = paste0(dataDir, 'scRNA-seq-regeneration.rds'))

DimPlot(aa, reduction = 'umap', raster = FALSE, group.by = 'annotation', label = TRUE, repel = TRUE) 


aa$celltypes = aa$annotation
aa$celltypes[grep('Cardiomyocytes', aa$celltypes)] = 'CMs'
aa$celltypes[grep('ECs', aa$celltypes)] = 'ECs'
aa$celltypes[grep('Endothelial', aa$celltypes)] = 'Endos'
aa$celltypes[grep('Endocardium', aa$celltypes)] = 'Endocardium'
aa$celltypes[grep('Fibroblasts', aa$celltypes)] = 'FBs'
aa$celltypes[grep('Macrophages', aa$celltypes)] = 'Macrophages'

p1 = DimPlot(aa, reduction = 'umap', raster = FALSE, group.by = 'celltypes', label = TRUE, repel = TRUE) 

p2 = FeaturePlot(aa, features = 'AXL', reduction = 'umap', raster = FALSE)

p1 / p2

ggsave(filename = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/Plots4manuscripts/revision_1/', 
                         'zebrafish_Li2025_Axl.pdf'), 
       width = 8, height = 12)


aa$condition = aa$time_points

p1 = DimPlot(aa, reduction = 'umap', raster = FALSE, group.by = 'condition', label = TRUE, repel = TRUE) 
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, 'zebrafish_Li2025_umap_condition_annotation.pdf'), 
       width = 20, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'zebrafish_Li2025_scRNAseq_regeneration_all.rds'))


## select only the CMs and plot AXLs
aa = subset(aa, cells = colnames(aa)[which(aa$celltypes == 'CMs')])

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.3)

aa$condition = aa$time_points

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE)
p3 = FeaturePlot(aa, features = c('AXL'))

(p1 + p2)/p3

ggsave(filename = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/Plots4manuscripts/revision_1/', 
                         'zebrafish_Li2025_CMs_Axl.pdf'), width = 12, height = 10)


##########################################
# an overview of zebrafish scRNA-seq data 
##########################################
aa = readRDS(file = paste0(RdataDir, 'zebrafish_Li2025_scRNAseq_regeneration_all.rds'))

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE)

p1 + p2

## discard the blood cells
aa = subset(aa, cells = colnames(aa)[which(aa$annotation != "Red blood cells")])

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.3)


