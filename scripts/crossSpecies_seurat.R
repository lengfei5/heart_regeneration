##########################################################################
##########################################################################
# Project: heart regeneration   
# Script purpose: test Seurat for cross-species maaping 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 12 10:27:53 2024
##########################################################################
##########################################################################
rm(list = ls())
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

outDir = paste0(resDir, '/seurat_out/')
system(paste0('mkdir -p ', outDir))

##########################################
# import ax, nm, mm data 
##########################################
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq.rds'))
nm = readRDS(file = paste0(RdataDir, 'nm_scRNAseq.rds'))
mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq.rds'))

aa = nm
aa$species = 'nm'
#aa$subtype = aa$FineID
#aa$celltype = aa$BroadID
#aa = subset(aa, cells = colnames(aa)[which(aa$celltype != 'glial')])
DefaultAssay(aa) = 'integrated'

aa@assays$integrated@counts = aa@assays$RNA@counts

aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'integrated', 
                dimreducs = c('mnn', 'umap'))

range(aa@assays$integrated@data)
xx = aa@assays$integrated@data
xx[which(xx<0)] = 0
aa@assays$integrated@data = xx
rm(xx)
range(aa@assays$integrated@data)
range(aa@assays$integrated@counts)

p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE) 
p3 = DimPlot(aa, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(aa, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)

mm$species = 'mm'

range(mm@assays$integrated@data)
xx = mm@assays$integrated@data
xx[which(xx<0)] = 0
mm@assays$integrated@data = xx
range(mm@assays$integrated@data)
range(mm@assays$integrated@counts)

rm(xx)

p1 = DimPlot(mm, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(mm, group.by = 'condition', label = TRUE, repel = TRUE) 
p3 = DimPlot(mm, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(mm, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)

ggs = intersect(rownames(aa), rownames(mm))
mm = subset(mm, features = ggs)
aa = subset(aa, features = ggs)

aa = merge(aa, y = mm, add.cell.ids = c("nm", "mm"), project = "heartReg")

range(aa@assays$integrated@data)

#saveRDS(aa, file = paste0(outDir, 'nm_mm_scRNAseq_merged_v1.rds'))

range(ax@assays$RNA@data)
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

new_ax <- CreateSeuratObject(counts=counts, assay = 'integrated', meta.data = metadata)
new_ax<- NormalizeData(new_ax, normalization.method = "LogNormalize", scale.factor = 10000)
new_ax <- FindVariableFeatures(new_ax, selection.method = "vst", nfeatures = 3000)

new_ax <- ScaleData(new_ax, features = rownames(new_ax))

range(new_ax@assays$integrated@data)

rm(list = c('counts', 'metadata', 'ax'))

aa = merge(aa, y = new_ax, add.cell.ids = c("m", "ax"), project = "heartReg")

rm(new_ax)
rm(mm)
rm(nm)

saveRDS(aa, file = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_forSeurat_v1.rds'))

##########################################
# test Seurat 
##########################################
source('functions_dataIntegration.R')
aa = readRDS(file = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_forSeurat_v1.rds'))

range(aa@assays$integrated@data)
aa = ScaleData(aa, features = rownames(aa))

method = 'Harmony'

if(method == 'Seurat_RPCA'){
  ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                           group.by = 'species', 
                                           nfeatures = 3000,
                                           merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                           redo.normalization.scaling = FALSE,
                                           correct.all = FALSE)
  
  p1 = DimPlot(ref.combined, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle('Seurat_RPCA')
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_RPCA")
  p3 = DimPlot(ref.combined, group.by = 'dataset', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_RPCA")
  
  p1 /(p2 +p3)
  
  ggsave(filename = paste0(outDir, '/cross_species_mapping_Seurat_RPCA.pdf'), 
         width = 16, height = 12)
  
  
  
}

if(method == 'Seurat_CCA'){
  ref.combined = IntegrateData_Seurat_CCA(aa, 
                                          group.by = 'species', 
                                          nfeatures = 3000,
                                          merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                          redo.normalization.scaling = FALSE,
                                          correct.all = FALSE)
  
  p1 = DimPlot(ref.combined, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle('Seurat_CCA')
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_CCA")
  p3 = DimPlot(ref.combined, group.by = 'dataset', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_CCA")
  
  p1 + p2
  
  ggsave(filename = paste0(outDir, '/cross_species_mapping_Seurat_CCA.pdf'), 
         width = 16, height = 6)
  
  saveRDS(ref.combined, file = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_SeuratCCA.rds'))
  
  
}


if(method == 'Harmony'){
  source('functions_dataIntegration.R')
  ref.combined = IntegrateData_runHarmony(aa, 
                                          group.by = 'species',
                                          nfeatures = 3000,
                                          dims.use = c(1:30),
                                          #merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                          #redo.normalization.scaling = FALSE,
                                          redo.normalization.hvg.scale.pca = FALSE,
                                          max.iter.harmony = 10
                                          #correct.all = FALSE
                                          )
  
  p1 = DimPlot(ref.combined, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle(method)
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle(method)
  p3 = DimPlot(ref.combined, group.by = 'dataset', label = TRUE, repel = TRUE) +
    ggtitle(method)
  
  p1 + p2
  
  ggsave(filename = paste0(outDir, '/cross_species_mapping_', method, '.pdf'), 
         width = 16, height = 6)
  
  #saveRDS(ref.combined, file = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_', method, '.rds'))
  
}

