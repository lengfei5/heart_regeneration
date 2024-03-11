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


########################################################
########################################################
# Section II: # add subtypes with the original subclusterd based on snRNA and scATAC
# original R objects from https://zenodo.org/records/7098004#.Y0P_LC0RoeY
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected.rds'))
aa$patient_region_id = droplevels(aa$patient_region_id)

aa$cell.id = sapply(colnames(aa), function(x) {unlist(strsplit(x, '_'))[1]})
aa$cell.id = paste0(aa$patient_region_id, '_', aa$cell.id)

# ax = readRDS(file = paste0())

DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

annot_file = list.files(path = '../published_dataset/human/Kuppe_et_al_2022/processed_data_Robj/subtype_annot',
                        pattern = '*.Rds', full.names = TRUE)

aa$annotation = NA

##########################################
# subtypes of CMs
##########################################
for(n in 1:length(annot_file))
{
  # n = 1
  celltype = gsub('_snRNA_snATAC.Rds', '', basename(annot_file[n]))
  cat(n, '--', basename(annot_file[n]), '\n')
  xx = readRDS(annot_file[n])
  cat(nrow(xx), 'cells \n')
  print(table(xx$annotation))
  
  xx$cell.id = sapply(colnames(xx), function(x) {unlist(strsplit(x, '#'))[2]})
  xx$cell.id = paste0(xx$patient_region_id, '_', xx$cell.id)
  
  xx = DietSeurat(xx, data = TRUE, assays = 'RNA')
  
  subs = subset(aa, cells = colnames(aa)[which(aa$cell_type_original == celltype)])
  DimPlot(subs, group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)
    
  #mm = match(xx$cell.id, aa$cell.id)
  #cat(length(which(!is.na(mm))), 'cell matches \n')
  subs <- FindVariableFeatures(aa, selection.method = 'vst', nfeatures = 2000)
  #subs <- ScaleData(subs)
  #subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = FALSE)
  
  #ElbowPlot(subs, ndims = 50)
  subs <- FindNeighbors(subs, dims = 1:20, reduction = 'harmony')
  subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  #Idents(subs) = subs$condition
  subs <- RunUMAP(subs, dims = 1:20, reduction = 'harmony', n.neighbors = 30, min.dist = 0.3)
  
  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/umap_timepoints_clusters_regress.nCounts.pdf'), width = 16, height = 6)
  
  
    
}


