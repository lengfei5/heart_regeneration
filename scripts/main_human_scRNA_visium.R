##########################################################################
##########################################################################
# Project: heart regeneration project
# Script purpose: process and analyze the scRNA-seq data and visium data from Kuppe_et_al_2022
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Feb 22 11:34:57 2024
##########################################################################
##########################################################################

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

p1 = DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

p2 = DimPlot(aa, reduction = 'umap', group.by = 'final_cluster', raster=FALSE, label = TRUE, repel = TRUE)

p1 + p2 

ggsave(paste0(resDir, '/Kupper2022_Umap_clusters_cellType.original.pdf'), 
       width = 16, height = 8)

p1 = DimPlot(aa, reduction = 'umap', group.by = 'patient_region_id', raster=FALSE, label = TRUE, repel = TRUE)
p2 = DimPlot(aa, reduction = 'umap', group.by = 'major_labl', raster=FALSE, label = TRUE, repel = TRUE)
p1 + p2

ggsave(paste0(resDir, '/Kupper2022_Umap_patienceID_major_labels.pdf'), 
       width = 16, height = 8)

saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_all.rds'))

##########################################
# subset the snRNA-seq with the donor based on the metadata (41586_2022_5060_MOESM6_ESM)
# donors with days after infarction >30 were filtered
##########################################
Idents(aa) = aa$donor_id

aa = subset(aa, idents = c("P14", "P18", 'P20', 'P19', 'P4', 'P13', 'P11', 'P12', 'P5'), invert = TRUE)

DimPlot(aa, reduction = 'umap', group.by = 'donor_id', raster=FALSE, label = TRUE, repel = TRUE)

ggsave(paste0(resDir, '/Kupper2022_Umap_selectedDonor.pdf'), 
       width = 16, height = 12)


p1 = DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

p2 = DimPlot(aa, reduction = 'umap', group.by = 'donor_id', raster=FALSE, label = TRUE, repel = TRUE)

p1 + p2 

ggsave(paste0(resDir, '/Kupper2022_Umap_clusters_cellType.original_selectedDonor.pdf'), 
       width = 16, height = 8)


saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected.rds'))

##########################################
# add subtypes with the original subclusterd based on snRNA and scATAC
# original R objects from https://zenodo.org/records/7098004#.Y0P_LC0RoeY
##########################################









