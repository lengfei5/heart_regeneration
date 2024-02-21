rm(list = ls())

version.analysis = '_20231207'

resDir = paste0("../results/scRNAseq_human", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

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
library(DropletUtils)
library(edgeR)
library(future)
options(future.globals.maxSize = 160000 * 1024^2)
library(pryr) # monitor the memory usage

mem_used()

Normalization = 'lognormal' # ('lognormal or SCT')

dataDir = "../published_dataset/human/Kuppe_et_al_2022/processed_data_Robj/snRNA/"

########################################################
########################################################
# Section I : Import the processed seurat object of scRNA-seq from
# Kuppe et al., 2022
# object were downloaded in https://cellxgene.cziscience.com/collections/8191c283-0816-424b-9b61-c3e1d6258a77
########################################################
########################################################
aa = readRDS(file = paste0(dataDir, 'all-snRNA.rds'))

DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

DimPlot(aa, reduction = 'umap', group.by = 'final_cluster', raster=FALSE, label = TRUE, repel = TRUE)

DimPlot(aa, reduction = 'umap', group.by = 'donor_id', raster=FALSE, label = TRUE, repel = TRUE)

saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_all.rds'))


