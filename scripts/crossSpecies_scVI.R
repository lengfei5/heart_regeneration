##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: prepare the scRNA-seq data for each species 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov 11 10:55:32 2024
##########################################################################
##########################################################################
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
library(pryr) # monitor the memory usage

source('functions_scRNAseq.R')
source('functions_Visium.R')

options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

version.analysis = '_20240220'
resDir = paste0("../results/cross_species", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

outDir = paste0(resDir, '/scvi_out/')
system(paste0('mkdir -p ', outDir))


##########################################
# import the neonatal mouse data before batch correction 
##########################################
aa = readRDS(file = paste0('../results/scRNAseq_neonatalMouse_20231108/Rdata/', 
                           'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_regress.nUMI.rds'))

aa$dataset = factor(aa$dataset, levels = c('Wang2020', 'Cui2020'))
aa$batch = paste0(aa$dataset, '_', aa$condition)
aa$batch = factor(aa$batch)

p1 = DimPlot(aa, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE) 
p3 = DimPlot(aa, group.by = 'BroadID', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(aa, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)


aa$species = 'nm'
aa$subtype = aa$FineID
aa$celltype = aa$BroadID

aa = subset(aa, cells = colnames(aa)[which(aa$celltype != 'glial')])


##########################################
# import adult mice data and merge with neonatal mouse 
##########################################
aa = readRDS(file = paste0('../results/Rdata/',
                           'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
cms = readRDS(file =  paste0('../results/Rdata/', 
                             'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap_subtypes.rds'))

aa$dataset = 'Forte2020'
cms$dataset = 'Ren2020'
#aa$annot.ref = aa$my_annot
#cms$annot.ref = cms$CellType
features.common = intersect(rownames(aa), rownames(cms))
length(union(rownames(aa), rownames(cms)))

mm = merge(aa, y = cms, add.cell.ids = c("Forte2020", "Ren2020"), project = "adultHeart")

#mm = readRDS(file = paste0('../results/Rdata/', 
#                           'Seurat.obj_adultMiceHeart_Forte2020_Ren2020_refCombined_logNormalize_',
#                           'counts_v3.rds'))

mm$species = 'mm'
mm$condition = mm$timepoints

mm = FindVariableFeatures(mm, selection.method = "vst", nfeatures = 3000)
mm <- ScaleData(mm, verbose = FALSE)
mm <- RunPCA(mm, npcs = 30, verbose = FALSE)

ElbowPlot(mm, ndims = 30)

mm <- FindNeighbors(mm, reduction = "pca", dims = 1:20)
#mm <- FindClusters(mm, resolution = 0.2)
mm <- RunUMAP(mm, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 


p1 = DimPlot(mm, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(mm, group.by = 'condition', label = TRUE, repel = TRUE) 
p3 = DimPlot(mm, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(mm, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)


ggs = intersect(rownames(aa), rownames(mm))

mm = subset(mm, features = ggs)
aa = subset(aa, features = ggs)

aa = merge(aa, y = mm, add.cell.ids = c("nm", "mm"), project = "heartReg")

saveRDS(aa, file = paste0(RdataDir, 'nm_mm_scRNAseq_merged_v1.rds'))

##########################################
# import the axolotl data  
##########################################
aa = readRDS(file = paste0(RdataDir, 'nm_mm_scRNAseq_merged_v1.rds'))

aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA')
table(aa$batch)

aa$dataset[which(aa$dataset == 'Cui2020')] = 'nm_Cui2020'
aa$dataset[which(aa$dataset == 'Wang2020')] = 'nm_Wang2020'
aa$dataset[which(aa$dataset == 'Forte2020')] = 'mm_Forte2020'
aa$dataset[which(aa$dataset == 'Ren2020')] = 'mm_Ren2020'

aa$batch[which(aa$dataset == 'mm_Forte2020')] = 'mm_Forte2020'
aa$batch[which(aa$dataset == 'mm_Ren2020')] = 'mm_Ren2020'

ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq.rds'))

ax$species = 'ax'
ax$dataset = 'ax_Bassat2024'
ax$batch = ax$dataset

ax$subtype = ax$celltypes

ax$celltype = ax$celltypes

ax$celltype[grep('FB_', ax$celltype)] = 'FB'
ax$celltype[grep('B_', ax$celltype)] = 'B'
ax$celltype[grep('CM_', ax$celltype)] = 'CM'
ax$celltype[grep('EC_', ax$celltype)] = 'EC'
ax$celltype[grep('Mo.Macs_', ax$celltype)] = 'Mo.Macs'

ax = subset(ax, cells = colnames(ax)[which(ax$celltype != 'Neuronal' & ax$celltype != 'RBC' &
                                             ax$celltype != 'Proliferating_RBC')])

ax = subset(ax, cells = colnames(ax)[which(ax$celltype != 'Neu_IL1R1')])

## define the one-on-one ortholog between axololt and mice
an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

jj = which(!is.na(match(an_orthologs$query, rownames(aa))))
an_orthologs = an_orthologs[jj, ]

counts = table(an_orthologs$query)
gg_uniq = names(counts)[which(counts == 1)]
jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

an_orthologs = an_orthologs[jj2, ]

ax = subset(ax, features = an_orthologs$ref)

aa = subset(aa, features = an_orthologs$query)


counts = ax@assays$RNA@counts
metadata = ax@meta.data
counts = counts[match(an_orthologs$ref, rownames(counts)), ]
rownames(counts) = an_orthologs$query

new_ax <- CreateSeuratObject(counts=counts, assay = 'RNA', meta.data = metadata)

rm(list = c('counts', 'metadata'))

aa = merge(aa, y = new_ax, add.cell.ids = c("m", "ax"), project = "heartReg")

rm(new_ax)
#aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA')
#jj = which(is.na(aa$batch))
#aa$batch[jj] = aa$dataset[]

saveRDS(aa, file = paste0(RdataDir, 'nm_mm_ax_scRNAseq_merged_v2.rds'))

library(SeuratDisk)
VariableFeatures(aa)
#Idents(aa) = droplevels(aa$condition)

saveFile = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_v1.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)

##########################################
# correct cell labels for scANVI
##########################################
aa = readRDS(paste0(RdataDir, 'nm_mm_ax_scRNAseq_merged_v2.rds'))

aa$celltype[which(aa$celltype == 'T_cells')] = 'T'
aa$celltype[which(aa$celltype == 'prolife.Mphage')] = 'MP'
aa$celltype[which(aa$celltype == 'MHCII.Mphage')] = 'MP'

aa$celltype[which(is.na(aa$celltype))] = 'Unknown'

aa$celltype[which(aa$subtype == 'CM_Cav3.1')] = 'Unknown'

library(SeuratDisk)
VariableFeatures(aa)
#Idents(aa) = droplevels(aa$condition)

saveFile = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_celllabels_v2.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)





