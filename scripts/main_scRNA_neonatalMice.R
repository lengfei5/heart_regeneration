##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: process and analyze the neonatal mice scRNA-seq/snRNA-seq
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Nov  8 13:29:36 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_20231108'

resDir = paste0("../results/scRNAseq_neonatalMouse", version.analysis)
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
library(pryr) # monitor the memory usage
mem_used()

Normalization = 'lognormal' # ('lognormal or SCT')

########################################################
########################################################
# Section I : First process the neonatal cardiomyocyte from 
# Cui et al., Dev Cell 2020
########################################################
########################################################
Process.Cardiomyocyte.Ren.2020 = FALSE
if(Process.Cardiomyocyte.Ren.2020){
  dataDir = '../published_dataset/neonatal_mice/Cui_2020_DevCell/'
  counts = fread(file = paste0(dataDir, 'GSE130699_P1_8_AllCM.txt'), nThread = 6)
  rownames(counts) = counts$V1
  
  #metadata = read.delim(file = paste0(dataDir, '/GSE120064_TAC_clean_cell_info_summary.txt'), 
  #                      sep = '\t', header = TRUE)
  
  # align the cellID in metadata and in the count table
  mm = match(metadata$CellID, colnames(counts))
  counts = counts[, c(1, mm), with=FALSE]
  
  # keep only the week 0 data
  jj = which(metadata$condition == '0w' | metadata$condition == '2w')
  
  metadata = metadata[jj, ]
  counts = counts[, c(1, jj+1), with = FALSE]
  counts = as.data.frame(counts)
  rownames(counts) = counts$V1
  counts = counts[, -1]
  
  # remove genes with 0 umi counts
  ss = apply(as.matrix(counts), 1, sum)
  counts = counts[which(ss>0), ]
  
  metadata = data.frame(metadata)
  rownames(metadata) = metadata$CellID
  
  save(metadata, counts, file = paste0(RdataDir, 'metadata_counts_week0.week2.Rdata'))
  
}

##########################################
# make Seurat object with metadata and counts  
##########################################
load(file = paste0(RdataDir, 'metadata_counts_week0.week2.Rdata'))

myo <- CreateSeuratObject(counts = counts, project = "adult", min.cells = 50, min.features = 500)
myo = AddMetaData(myo, metadata, col.name = NULL) 

myo[["percent.mt"]] <- PercentageFeatureSet(myo, pattern = "^mt-")

plot1 <- FeatureScatter(myo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(myo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

aa = subset(myo, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 5*10^5)
rm(myo)
rm(counts)

if(SCT.normalization){
  aa <- SCTransform(aa, assay = "RNA", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE, 
                    min_cells=5) 
  
}else{
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  plot1 <- VariableFeaturePlot(aa)
  
  top10 <- head(VariableFeatures(aa), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  
}

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:10)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

### KEEP the umap configuration for only week0: nfeatures = 2000, dims = 1:20, n.neighbors = 30 and min.dist = 0.05 (sctransform)
### umap configuraiton for week0 and week2: nfeatures = 3000, dims = 1:20, n.neighbors = 30/50, min.dist = 0.1 (sctransform)
### umap configuration for week0 and week2: nfeatures = 3000, dims = 1:20, n.neighbors = 20, min.dist = 0.05 (seurat.norm)

if(Normalization == 'SCT'){
  aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
}else{
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
  #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
}
