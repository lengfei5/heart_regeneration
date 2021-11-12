##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: prepare the single cell/nucleus RNA-seq data for the cell type reference
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Nov 11 14:14:12 2021
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_scRNAseq_reference_20211111'

resDir = paste0("../results/scRNAseq_mouse", version.analysis)
RdataDir = paste0('../results/Rdata/')

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

########################################################
########################################################
# Section : First process the adult cardiomyocyte from Ren et al., 2020
# 
########################################################
########################################################
dataDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/published_dataset/Ren_2020'
metadata = read.delim(file = paste0(dataDir, '/GSE120064_TAC_clean_cell_info_summary.txt'), sep = '\t', header = TRUE)
counts = fread(file = paste0(dataDir, '/GSE120064_TAC_raw_umi_matrix.csv'), header = TRUE, nThread = 6)
#rownames(counts) = counts$V1

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

save(metadata, counts, file = paste0(RdataDir, 'metadata_counts_week.0.2.Rdata'))

##########################################
# make Seurat object with metadata and counts  
##########################################
myo <- CreateSeuratObject(counts = counts, project = "adult", min.cells = 50, min.features = 500)
myo = AddMetaData(myo, metadata, col.name = NULL) 

myo[["percent.mt"]] <- PercentageFeatureSet(myo, pattern = "^mt-")

plot1 <- FeatureScatter(myo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(myo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

myo = subset(myo, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 5*10^5)

aa = myo

aa <- SCTransform(aa, assay = "RNA", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE, 
                  min_cells=5) 

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:10)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

####### KEEP this umap configuration: nfeatures = 2000, dims = 1:20, n.neighbors = 30 and min.dist = 0.05 for week0 dataset
#######
aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.1) 

p1 = DimPlot(aa, reduction = 'umap', group.by = 'seurat_clusters')
p2 = DimPlot(aa, reduction = "umap", group.by = c("CellType"))
p1 + p2 


p1 = DimPlot(aa, reduction = 'umap', group.by = 'condition')
p2 = DimPlot(aa, reduction = 'umap', group.by = 'sample')
p1 + p2

gg.examples =  c('Nppa', 'Nppb', 'Myh6', 'Tnnc1', 'Tnni3', 
                 'Tnnt2', 'Actn2', 'Gata4', 'Nkx2-5', 'Gja1', 
                 'Myl2', 'Tpm1', 'Ryr2', 'Atp2a2', 'Acta1')

p3 = FeaturePlot(aa, reduction = 'umap', features = gg.examples)


p2 + p3 
p1 + p3

p1 + p2 + ggsave(paste0(resDir, '/Umap_newClusters_vs_cellType.original_v2.pdf'), 
                 width = 12, height = 8)

#DimPlot(aa, reduction = "umap", group.by = c("SubCluster"))

FeaturePlot(aa, reduction = 'umap', features = c('Csf1r', 'Cd163'))
FeaturePlot(aa, reduction = 'umap', features = c('S100a3', 'S100a9'))

##########################################
# heatmap of marker genes
##########################################
# remove the mitochrio marker genes
jj = grep('^mt-', rownames(aa)) 
aa = aa[-jj, ]


aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
aa.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(aa, features = top10$gene) + NoLegend() + ggsave(paste0(resDir, '/heatmap_markerGenes_rmMt.pdf'), 
                                                           width = 12, height = 16)

