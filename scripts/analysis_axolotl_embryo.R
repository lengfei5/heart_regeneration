##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: analyze the axolotl cardiac_lineage data from Sebastian
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep 14 13:33:26 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_20260914'

resDir = paste0("../results/axoltol_embryo_cardiac_lineage", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/Collaborations/Jingkui-Elad/From_Sebestian_embryo'

source('functions_scRNAseq.R')
source('functions_Visium.R')
species = 'axloltl_scRNAseq'

require(Seurat)
#require(sctransform)3
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 80000 * 1024^2)

mem_used()


########################################################
########################################################
# Section I: double check the QC and overview of the data
# 
########################################################
########################################################

##########################################
# this option did not work
##########################################
library(Seurat)
library(SeuratDisk)
# Convert h5ad -> h5seurat first
Convert(paste0(dataDir, "/adata_cardiac_lineage_subset.h5ad"), dest = "h5seurat", 
        overwrite = TRUE)

# Load the h5seurat file
seurat_obj <- LoadH5Seurat("data.h5seurat")


##########################################
# try sceasy option, didn't work either
##########################################
library(sceasy)
library(reticulate)

use_condaenv('scvi-env')

sceasy::convertFormat(paste0(dataDir, "/adata_cardiac_lineage_subset.h5ad"), 
                      from = "anndata", to = "seurat",
                      outFile = "data.rds")



##########################################
# import the count table, metadata and also dimension reduction 
##########################################
counts = read.csv(paste0(RdataDir, '/cardiac_lineage_subset_adata_layer_counts.csv'), 
                  header = TRUE, sep = ',', row.names = c(1))
metadata = read.csv(paste0(RdataDir, '/cardiac_lineage_subset_adata_metadata.csv'), 
                    header = TRUE, sep = ',', row.names = c(1))
pca = read.csv(paste0(RdataDir, 'cardiac_lineage_subset_adata_rd_pca.csv'), 
               header = TRUE, sep = ',', row.names = c(1))
umap = read.csv(paste0(RdataDir, 'cardiac_lineage_subset_adata_rd_umap.csv'), 
                header = TRUE, sep = ',', row.names = c(1))

aa <- CreateSeuratObject(counts = t(as.matrix(counts)), project = "cardiacEmbryo", assay = "RNA",
                         min.cells = 3, min.features = 50, meta.data = metadata)

colnames(umap) <- paste0("UMAP_", 1:ncol(umap))
aa[['umap_orig']] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAPorig_", 
                                          assay = DefaultAssay(aa))

colnames(pca) = paste0("PCA_", 1:ncol(pca))
aa[['pca_orig']] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PCAorig_", 
                                         assay = DefaultAssay(aa))
#rm(list(c(counts, metadata, pca, umap)))
rm(counts); rm(pca); rm(umap)

saveRDS(aa, file = paste0(RdataDir, 'axolotlEmbryo_cardiac_lineage_subset.rds'))


##########################################
# processing the counts
##########################################
aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 30)
aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(aa, label = TRUE, reduction = 'umap_orig', repel = TRUE, group.by = 'stage')

p2 = DimPlot(aa, label = TRUE, reduction = 'umap_orig', repel = TRUE, group.by = 'cell_type')

p1

ggsave(filename = paste0(resDir, '/axolotlEmbryo_cardiakLineage_umap_stage.pdf'), 
       width = 12, height = 8)

p2

ggsave(filename = paste0(resDir, '/axolotlEmbryo_cardiakLineage_umap_celltypes.pdf'), 
       width = 12, height = 8)


mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4", "ATP8", "MT-CO1", "COI")
mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
mtgenes = rownames(aa)[!is.na(match(rownames(aa), mtgenes))]

xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "RNA", features = mtgenes)
aa[['percent.mt']] = xx$percent.mt

rm(xx)


FeaturePlot(aa, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'))

ggsave(filename = paste0(resDir, '/axolotlEmbryo_cardiakLineage_umap_QCs.pdf'), 
       width = 12, height = 16)

VlnPlot(aa ,features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'stage')

ggsave(filename = paste0(resDir, '/axolotlEmbryo_cardiakLineage_vlnplot_QCs.pdf'), 
       width = 20, height = 8)
